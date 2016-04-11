/*
 * MatrixFactorization.cpp
 *
 *  Created on: Jun 29, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "MatrixFactorization.hpp"

#include <boost/assign.hpp>
#include <fstream>
#include <string>

#include <boost/assign.hpp>

#include "Log/Logging.hpp"
#include "MatrixFactorizer/Database/IterationInformation.hpp"
#include "MatrixFactorizer/Helpers/detail.hpp"
#include "MatrixFactorizer/Options/MatrixFactorizerOptions.hpp"
#include "MatrixFactorizer/Work/Master.hpp"
#include "MatrixFactorizer/ScalingAmbiguityRemover/MatrixProductEqualizer.hpp"
#include "MatrixFactorizer/ScalingAmbiguityRemover/MatrixProductRenormalizer.hpp"
#include "MatrixFactorizer/ScalingAmbiguityRemover/ScalingAmbiguityMaintainer.hpp"
#include "MatrixFactorizer/Solvers/InRangeSolver.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/VectorSetter.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriteriaFactory.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingArguments.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraintsFactory.hpp"

using namespace boost::assign;

MatrixFactorization::MatrixFactorization(
		const MatrixFactorizerOptions &_opts,
		IterationInformation &_info
#ifdef MPI_FOUND
		,boost::mpi::communicator &_world
#endif
		) :
		pixel_opts(_opts),
		spectral_opts(_opts),
		opts(_opts),
		info(_info)
#ifdef MPI_FOUND
		, world(_world)
#endif
{
	// an operator between two lp spaces always factors through a
	// Hilbert space. Hence, we use an l_2 space between the two
	// matrix factors
	const_cast<CommandLineOptions &>(pixel_opts).type_spacey =
			"lp";
	const_cast<CommandLineOptions &>(pixel_opts).py = 2.;
	const_cast<CommandLineOptions &>(pixel_opts).powery = 2.;
	const_cast<CommandLineOptions &>(spectral_opts).type_spacex =
			"lp";
	const_cast<CommandLineOptions &>(spectral_opts).px = 2.;
	const_cast<CommandLineOptions &>(spectral_opts).powerx = 2.;
}

void MatrixFactorization::operator()(
		const Eigen::MatrixXd &_data,
		int &_returnstatus
		)
{
#ifdef MPI_FOUND
	Master master(world, opts.overall_keys);

	// So far, Slaves are present and expect initial go (or not full_terminate
	// signal). Hence, if something has gone wrong, then we need at least to
	// tell them here before the actual solver loop
	if (_returnstatus != 0) {
		master.sendTerminate();
	} else
#endif
	{
		/// start timing
		boost::chrono::high_resolution_clock::time_point timing_start =
				boost::chrono::high_resolution_clock::now();

		unsigned int loop_nr = 0;

		// calculate an initial residual prior to the starting factors
		double residual = _data.norm();
		BOOST_LOG_TRIVIAL(info)
			<< "#" << loop_nr << " starting residual is " << residual;
		info.replace(IterationInformation::LoopTable, "residual", residual);
		info.addTuple(IterationInformation::LoopTable);

		/// parse in factors
		bool parse_status = true;
		if (opts.DoParseFactors) {
			const int result_parse_spectralmatrix = detail::parseFactorFile(
					opts.solution_factor_one_file.string(), spectral_matrix, "first factor");
			parse_status &= (result_parse_spectralmatrix == 0);
			const int result_parse_pixelmatrix = detail::parseFactorFile(
					opts.solution_factor_two_file.string(), pixel_matrix, "second factor");
			parse_status &= (result_parse_pixelmatrix == 0);

			// check whether parsing was successful, break otherwise
			if (!parse_status) {
				BOOST_LOG_TRIVIAL(error)
						<< "Parsing of either or both factors failed, aborting.";
				_returnstatus = 255;
#ifdef MPI_FOUND
				master.sendTerminate();
#endif
				return;
			}
		} else
			parse_status = false;

		/// construct solution starting points
		if (!parse_status) {
			BOOST_LOG_TRIVIAL(info)
					<< "Setting spectral matrix K to random starting values.";
			spectral_matrix = Eigen::MatrixXd(_data.rows(), opts.sparse_dim);
			detail::constructRandomMatrix(spectral_matrix);
		} else {
			if ((spectral_matrix.rows() == _data.rows())
					&& (spectral_matrix.cols() == opts.sparse_dim)) {
				BOOST_LOG_TRIVIAL(info)
						<< "Using spectral matrix K parsed from file.";
			} else {
				BOOST_LOG_TRIVIAL(error)
						<< "Parsed spectral matrix has dimensions "
						<< spectral_matrix.rows() << "," << spectral_matrix.cols()
						<< ", while expecting " << _data.rows() << "," << opts.sparse_dim;
				_returnstatus = 255;
#ifdef MPI_FOUND
				master.sendTerminate();
#endif
				return;
			}
		}
		if (!parse_status) {
			BOOST_LOG_TRIVIAL(info)
					<< "Setting pixel matrix X to zero.";
			pixel_matrix = Eigen::MatrixXd(opts.sparse_dim, _data.cols());
			detail::constructZeroMatrix(
					pixel_matrix);
		} else {
			if ((pixel_matrix.rows() == opts.sparse_dim)
					&& (pixel_matrix.cols() == _data.cols())) {
				BOOST_LOG_TRIVIAL(info)
						<< "Using pixel matrix X parsed from file.";
			} else {
				BOOST_LOG_TRIVIAL(error)
						<< "Parsed pixel matrix has dimensions "
						<< pixel_matrix.rows() << "," << pixel_matrix.cols()
						<< ", while expecting " << opts.sparse_dim << "," <<  _data.cols();
				_returnstatus = 255;
#ifdef MPI_FOUND
				master.sendTerminate();
#endif
				return;
			}
		}

		// create outer ("loop") stopping criteria
		StoppingArguments stopping_args;
		stopping_args.setTolerance(opts.residual_threshold);
		stopping_args.setMaxIterations(opts.max_loops);
		std::string stopping_criteria(
				"MaxIterationCount || RelativeChangeResiduum || Residuum");
		StoppingCriteriaFactory stop_factory;
		StoppingCriterion::ptr_t stopping_criterion =
				stop_factory.create(stopping_criteria, stopping_args);

		// create auxiliary constraints
		AuxiliaryConstraintsFactory constraint_factory;
		AuxiliaryConstraints::ptr_t auxiliary_constraints =
				constraint_factory.create(opts.auxiliary_constraints);

		/// create feasible solution starting points
		if (auxiliary_constraints) {
			if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
				BOOST_LOG_TRIVIAL(trace)
						<< "Spectral matrix K^t before constraining is \n" << spectral_matrix.transpose();
			} else {
				BOOST_LOG_TRIVIAL(info)
						<< "Spectral matrix K^t before constraining is \n" << spectral_matrix.transpose();
			}

			NormedSpaceFactory::args_t args_SpaceL2;
			args_SpaceL2 +=
					boost::any(2.),
					boost::any(2.);
			NormedSpace_ptr_t L2 =
					NormedSpaceFactory::create(
							opts.sparse_dim, "lp", args_SpaceL2);
			spectral_matrix.transposeInPlace();
			for (int i=0;i<spectral_matrix.cols();++i) {
				SpaceElement_ptr_t tempvector = ElementCreator::create(
						L2,
						spectral_matrix.col(i));
//				BOOST_LOG_TRIVIAL(info)
//							<< "tempvector of col " << i << " before constraining is " << *tempvector;
				(*auxiliary_constraints)(tempvector);
//				BOOST_LOG_TRIVIAL(info)
//							<< "tempvector of col " << i << " after constraining is " << *tempvector;
				const Eigen::Ref<Eigen::MatrixXd::ColXpr> col = spectral_matrix.col(i);
				VectorSetter::set(col, tempvector);
			}
			spectral_matrix.transposeInPlace();

			if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
				BOOST_LOG_TRIVIAL(trace)
							<< "Spectral matrix K^t after constraining is \n" << spectral_matrix.transpose();
			} else {
				BOOST_LOG_TRIVIAL(info)
							<< "Spectral matrix K^t after constraining is \n" << spectral_matrix.transpose();
			}
		}

		// evaluate stop condition
		residual = detail::calculateResidual(_data, spectral_matrix, pixel_matrix);
		bool stop_condition = (*stopping_criterion)(
				boost::chrono::duration<double>(0),
				loop_nr,
				residual,
				1.);

		/// iterate over the two factors
		while (!stop_condition) {
			// update loop count
			++loop_nr;

			if (opts.fix_factor != 1) {
				info.replace(IterationInformation::LoopTable, "loop_nr", 2*(int)(loop_nr-1)+1);
				BOOST_LOG_TRIVIAL(debug)
					<< "======================== #" << loop_nr << "/1 ==================";

		//		renormalizeMatrixByTrace(spectral_matrix);

				if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
					BOOST_LOG_TRIVIAL(trace)
							<< "Current spectral matrix K^t is\n" << spectral_matrix.transpose();
				} else {
					BOOST_LOG_TRIVIAL(info)
							<< "Current spectral matrix K^t is\n" << spectral_matrix.transpose();
				}

#ifdef MPI_FOUND
				if (world.size() == 1) {
#endif
					if (!stop_condition) {
						InRangeSolver solver(
								spectral_opts,
								opts.overall_keys,
								opts.projection_delta);
						for (unsigned int dim = 0;
								(!stop_condition) && (dim < _data.cols());
								++dim) {
							Eigen::VectorXd solution;
							const bool solver_ok =
									solver(spectral_matrix,
											_data.col(dim),
											pixel_matrix.col(dim),
											solution,
											dim,
											auxiliary_constraints
											);
							if (solver_ok)
								pixel_matrix.col(dim) = solution;
							else
								BOOST_LOG_TRIVIAL(error)
									<< "The minimizer for InRangeSolver on spectralmatrix could not finish.";
							stop_condition |= !solver_ok;
						}
						// place accumulated values in loop table
						solver.insertAccumulatedProjectorValues(
								info.getLoopTable(), "_projection");
						solver.insertAccumulatedSolverValues(
								info.getLoopTable(), "_minimization");
					}
#ifdef MPI_FOUND
				} else {
					if (!stop_condition) {
						const bool solver_ok =
								master.solve(
										spectral_opts,
										spectral_matrix,
										_data,
										pixel_matrix,
										opts.auxiliary_constraints
										);
						// place accumulated values in loop table
						master.insertAccumulatedProjectorValues(
								info.getLoopTable(), "_projection");
						master.resetAccumulatedProjectorValues();
						master.insertAccumulatedSolverValues(
								info.getLoopTable(), "_minimization");
						master.resetAccumulatedSolverValues();
						if (!solver_ok)
							BOOST_LOG_TRIVIAL(error)
								<< "The minimizer for InRangeSolver on spectralmatrix could not finish.";
						stop_condition |= !solver_ok;
					}
				}
#endif

				if ((pixel_matrix.innerSize() > 10) || (pixel_matrix.outerSize() > 10)) {
					BOOST_LOG_TRIVIAL(trace)
							<< "Resulting pixel matrix X is\n" << pixel_matrix;
				} else {
					BOOST_LOG_TRIVIAL(info)
							<< "Resulting pixel matrix X is\n" << pixel_matrix;
				}

				// check criterion
				{
					residual =
							detail::calculateResidual(_data, spectral_matrix, pixel_matrix);
					info.replace(IterationInformation::LoopTable, "residual", residual);
					BOOST_LOG_TRIVIAL(info)
						<< "#" << loop_nr << " 1/2, residual is " << residual;
				}

				// submit loop tuple
				info.addTuple(IterationInformation::LoopTable);
			}

			if (opts.fix_factor != 2) {
				BOOST_LOG_TRIVIAL(debug)
					<< "======================== #" << loop_nr << "/2 ==================";

				info.replace(IterationInformation::LoopTable, "loop_nr", 2*(int)(loop_nr));

		//		renormalizeMatrixByTrace(pixel_matrix);

				if ((pixel_matrix.innerSize() > 10) || (pixel_matrix.outerSize() > 10)) {
					BOOST_LOG_TRIVIAL(trace)
							<< "Current pixel matrix X is\n" << pixel_matrix;
				} else {
					BOOST_LOG_TRIVIAL(info)
							<< "Current pixel matrix X is\n" << pixel_matrix;
				}

				// must transpose in place, as spectral_matrix.transpose() is const
				spectral_matrix.transposeInPlace();
#ifdef MPI_FOUND
				if (world.size() == 1) {
# endif
					if (!stop_condition) {
						InRangeSolver solver(
								pixel_opts,
								opts.overall_keys,
								opts.projection_delta);
						for (unsigned int dim = 0;
								(!stop_condition) && (dim < _data.rows());
								++dim) {
							Eigen::VectorXd solution;
							const bool solver_ok =
									solver(pixel_matrix.transpose(),
											_data.row(dim).transpose(),
											spectral_matrix.col(dim),
											solution,
											dim,
											auxiliary_constraints
											);
							if (solver_ok)
								spectral_matrix.col(dim) = solution;
							else
								BOOST_LOG_TRIVIAL(error)
									<< "The minimizer for InRangeSolver on pixelmatrix could not finish.";
							stop_condition |= !solver_ok;
						}
						// place accumulated values in loop table
						solver.insertAccumulatedProjectorValues(
								info.getLoopTable(), "_projection");
						solver.insertAccumulatedSolverValues(
								info.getLoopTable(), "_minimization");
					}
#ifdef MPI_FOUND
				} else {
					if (!stop_condition) {
						const bool solver_ok =
								master.solve(
										pixel_opts,
										pixel_matrix.transpose(),
										_data.transpose(),
										spectral_matrix,
										opts.auxiliary_constraints
										);
						// place accumulated values in loop table
						master.insertAccumulatedProjectorValues(
								info.getLoopTable(), "_projection");
						master.resetAccumulatedProjectorValues();
						master.insertAccumulatedSolverValues(
								info.getLoopTable(), "_minimization");
						master.resetAccumulatedSolverValues();
						if (!solver_ok)
							BOOST_LOG_TRIVIAL(error)
								<< "The minimizer for InRangeSolver on spectralmatrix could not finish.";
						stop_condition |= !solver_ok;
					}
				}
#endif
				spectral_matrix.transposeInPlace();

				// check criterion
				{
					residual =
							detail::calculateResidual(_data, spectral_matrix, pixel_matrix);
					BOOST_LOG_TRIVIAL(info)
						<< "#" << loop_nr << " 2/2, residual is " << residual;
				}
			}

			// remove ambiguity
			if (opts.fix_factor == 0) {
				MatrixProductEqualizer equalizer;
				const double scaling_change =
						equalizer(spectral_matrix, pixel_matrix);
				info.replace(IterationInformation::LoopTable, "scaling_change", scaling_change);
				BOOST_LOG_TRIVIAL(info)
					<< "Scaling factor is " << scaling_change;
			}

			if (opts.fix_factor != 2) {
				if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
					BOOST_LOG_TRIVIAL(trace)
							<< "Resulting spectral K^t matrix is\n" << spectral_matrix.transpose();
				} else {
					BOOST_LOG_TRIVIAL(info)
							<< "Resulting spectral K^t matrix is\n" << spectral_matrix.transpose();
				}

				// check criterion
				{
					residual =
							detail::calculateResidual(_data, spectral_matrix, pixel_matrix);
					BOOST_LOG_TRIVIAL(info)
						<< "#" << loop_nr << ", residual is " << residual;
					info.replace(IterationInformation::LoopTable, "residual", residual);
				}

				// submit loop tuple
				info.addTuple(IterationInformation::LoopTable);
			}

			stop_condition |= (*stopping_criterion)(
					boost::chrono::duration<double>(0),
					loop_nr,
					residual,
					1.);
		}
		if (loop_nr > opts.max_loops)
			BOOST_LOG_TRIVIAL(error)
				<< "Maximum number of loops " << opts.max_loops
				<< " exceeded, stopping iteration.";
		else
			BOOST_LOG_TRIVIAL(info)
				<< "Loop iteration performed " << loop_nr
				<< " times.";

		/// end timing
		boost::chrono::high_resolution_clock::time_point timing_end =
				boost::chrono::high_resolution_clock::now();
		BOOST_LOG_TRIVIAL(info) << "The operation took "
				<< boost::chrono::duration<double>(timing_end - timing_start).count()
				<< " seconds.";

		/// set last iteration information entry
		info.replace(IterationInformation::OverallTable, "loops", (int)loop_nr);
		info.replace(IterationInformation::OverallTable, "residual", residual);
		info.replace(IterationInformation::OverallTable, "runtime",
				boost::chrono::duration<double>(timing_end - timing_start).count());
		info.addTuple(IterationInformation::OverallTable);

#ifdef MPI_FOUND
		master.sendTerminate();
#endif

		// renormalize both factors
		MatrixProductRenormalizer renormalizer;
		renormalizer(spectral_matrix, pixel_matrix);

	}
}
