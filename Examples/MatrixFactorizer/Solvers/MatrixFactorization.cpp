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
#include <numeric>
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
#endif /* MPI_FOUND */
		) :
		pixel_opts(_opts),
		spectral_opts(_opts),
		info(_info)
#ifdef MPI_FOUND
		, world(_world)
#endif /* MPI_FOUND */
{
	// an operator between two lp spaces always factors through a
	// Hilbert space. Hence, we use an l_2 space between the two
	// matrix factors
	const_cast<MatrixFactorizerOptions &>(pixel_opts).type_spacey =
			"lp";
	const_cast<MatrixFactorizerOptions &>(pixel_opts).py = 2.;
	const_cast<MatrixFactorizerOptions &>(pixel_opts).powery = 2.;
	const_cast<MatrixFactorizerOptions &>(spectral_opts).type_spacex =
			"lp";
	const_cast<MatrixFactorizerOptions &>(spectral_opts).px = 2.;
	const_cast<MatrixFactorizerOptions &>(spectral_opts).powerx = 2.;
}

#ifdef MPI_FOUND
static void solveOneLoop_MPI(
		const Eigen::MatrixXd &_data,
		const Eigen::MatrixXd &fixed_factor,
		Eigen::MatrixXd &variable_factor,
		Master &master,
		const MatrixFactorizerOptions &opts,
		bool &stop_condition,
		IterationInformation &info
		)
{
	const bool solver_ok =
			master.solve(
					opts,
					fixed_factor,
					_data,
					variable_factor
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
			<< "The minimizer for InRangeSolver could not finish.";
	stop_condition |= !solver_ok;
}

#endif /* MPI_FOUND */

#ifdef OPENMP_FOUND
static void solveOneLoop_OpenMP(
		const Eigen::MatrixXd &_data,
		const Eigen::MatrixXd &fixed_factor,
		Eigen::MatrixXd &variable_factor,
		const MatrixFactorizerOptions &opts,
		bool &stop_condition,
		const AuxiliaryConstraints::ptr_t &auxiliary_constraints,
		IterationInformation &info
		)
{
	std::vector<double> returnvalues(_data.cols(), false);
	if (!stop_condition) {
#pragma omp parallel shared(returnvalues)
		{
			InRangeSolver solver(
					opts,
					opts.overall_keys,
					opts.projection_delta);
#pragma omp for
			for (unsigned int dim = 0; dim < _data.cols(); ++dim) {
				Eigen::VectorXd solution;
				returnvalues[dim] =
					solver(fixed_factor,
							_data.col(dim),
							variable_factor.col(dim),
							solution,
							dim,
							auxiliary_constraints
							);
				if (returnvalues[dim]) {
					#pragma omp critical (solution)
					variable_factor.col(dim) = solution;
				} else {
					BOOST_LOG_TRIVIAL(error)
						<< "The minimizer for InRangeSolver could not finish on column "
						<< dim << ".";
				}
			}
			// place accumulated values in loop table
#pragma omp single nowait
			{
				solver.insertAccumulatedProjectorValues(
						info.getLoopTable(), "_projection");
				solver.insertAccumulatedSolverValues(
						info.getLoopTable(), "_minimization");
			}
		}
		stop_condition |=
				std::find(returnvalues.begin(), returnvalues.end(), false) != returnvalues.end();
	}
}
#else /* OPENMP_FOUND */
static void solveOneLoop_sequentially(
		const Eigen::MatrixXd &_data,
		const Eigen::MatrixXd &fixed_factor,
		Eigen::MatrixXd &variable_factor,
		const MatrixFactorizerOptions &opts,
		bool &stop_condition,
		const AuxiliaryConstraints::ptr_t &auxiliary_constraints,
		IterationInformation &info
		)
{
	InRangeSolver solver(
			opts,
			opts.overall_keys,
			opts.projection_delta);
	for (unsigned int dim = 0;
			(!stop_condition) && (dim < _data.cols());
			++dim) {
		Eigen::VectorXd solution;
		const bool solver_ok =
		solver(fixed_factor,
				_data.col(dim),
				variable_factor.col(dim),
				solution,
				dim,
				auxiliary_constraints
				);
		if (solver_ok)
			variable_factor.col(dim) = solution;
		else
			BOOST_LOG_TRIVIAL(error)
				<< "The minimizer for InRangeSolver could not finish on column "
				<< dim << ".";
		stop_condition |= !solver_ok;
	}
	// place accumulated values in loop table
	solver.insertAccumulatedProjectorValues(
			info.getLoopTable(), "_projection");
	solver.insertAccumulatedSolverValues(
			info.getLoopTable(), "_minimization");
}
#endif /* OPENMP_FOUND */

void MatrixFactorization::operator()(
		const Eigen::MatrixXd &_data,
		int &_returnstatus
		)
{
#ifdef MPI_FOUND
	Master master(world, spectral_opts.overall_keys);
#endif /* MPI_FOUND */
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
		if (spectral_opts.DoParseFactors) {
			const int result_parse_spectralmatrix = detail::parseFactorFile(
					spectral_opts.solution_factor_one_file.string(), spectral_matrix, "first factor");
			parse_status &= (result_parse_spectralmatrix == 0);
			const int result_parse_pixelmatrix = detail::parseFactorFile(
					spectral_opts.solution_factor_two_file.string(), pixel_matrix, "second factor");
			parse_status &= (result_parse_pixelmatrix == 0);

			// check whether parsing was successful, break otherwise
			if (!parse_status) {
				BOOST_LOG_TRIVIAL(error)
						<< "Parsing of either or both factors failed, aborting.";
				_returnstatus = 255;
#ifdef MPI_FOUND
				master.sendTerminate();
#endif /* MPI_FOUND */
				return;
			}
		} else
			parse_status = false;

		/// construct solution starting points
		if (!parse_status) {
			BOOST_LOG_TRIVIAL(info)
					<< "Setting spectral matrix K to random starting values.";
			spectral_matrix = Eigen::MatrixXd(_data.rows(), spectral_opts.sparse_dim);
			detail::constructRandomMatrix(spectral_matrix);
		} else {
			if ((spectral_matrix.rows() == _data.rows())
					&& (spectral_matrix.cols() == spectral_opts.sparse_dim)) {
				BOOST_LOG_TRIVIAL(info)
						<< "Using spectral matrix K parsed from file.";
			} else {
				BOOST_LOG_TRIVIAL(error)
						<< "Parsed spectral matrix has dimensions "
						<< spectral_matrix.rows() << "," << spectral_matrix.cols()
						<< ", while expecting " << _data.rows() << "," << spectral_opts.sparse_dim;
				_returnstatus = 255;
#ifdef MPI_FOUND
				master.sendTerminate();
#endif /* MPI_FOUND */
				return;
			}
		}
		if (!parse_status) {
			BOOST_LOG_TRIVIAL(info)
					<< "Setting pixel matrix X to zero.";
			pixel_matrix = Eigen::MatrixXd(spectral_opts.sparse_dim, _data.cols());
			detail::constructZeroMatrix(
					pixel_matrix);
		} else {
			if ((pixel_matrix.rows() == spectral_opts.sparse_dim)
					&& (pixel_matrix.cols() == _data.cols())) {
				BOOST_LOG_TRIVIAL(info)
						<< "Using pixel matrix X parsed from file.";
			} else {
				BOOST_LOG_TRIVIAL(error)
						<< "Parsed pixel matrix has dimensions "
						<< pixel_matrix.rows() << "," << pixel_matrix.cols()
						<< ", while expecting " << spectral_opts.sparse_dim << "," <<  _data.cols();
				_returnstatus = 255;
#ifdef MPI_FOUND
				master.sendTerminate();
#endif /* MPI_FOUND */
				return;
			}
		}

		// create outer ("loop") stopping criteria
		StoppingArguments stopping_args;
		stopping_args.setTolerance(spectral_opts.residual_threshold);
		stopping_args.setMaxIterations(spectral_opts.max_loops);
		std::string stopping_criteria(
				"MaxIterationCount || RelativeChangeResiduum || Residuum");
		StoppingCriteriaFactory stop_factory;
		StoppingCriterion::ptr_t stopping_criterion =
				stop_factory.create(stopping_criteria, stopping_args);

		// create auxiliary constraints
		AuxiliaryConstraintsFactory constraint_factory;
		AuxiliaryConstraints::ptr_t auxiliary_constraints =
				constraint_factory.create(spectral_opts.auxiliary_constraints);

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
							spectral_opts.sparse_dim, "lp", args_SpaceL2);
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

			if (spectral_opts.fix_factor != 1) {
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

				if (!stop_condition) {
#ifdef MPI_FOUND
					if (world.size() != 1) {
						solveOneLoop_MPI(
								_data,
								spectral_matrix,
								pixel_matrix,
								master,
								spectral_opts,
								stop_condition,
								info);
					} else {
#else /* MPI_FOUND */
					{
#endif /* MPI_FOUND */
#ifdef OPENMP_FOUND
						solveOneLoop_OpenMP(
#else /* OPENMP_FOUND */
						solveOneLoop_sequentially(
#endif /* OPENMP_FOUND */
								_data,
								spectral_matrix,
								pixel_matrix,
								spectral_opts,
								stop_condition,
								auxiliary_constraints,
								info);
					}
				}

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

			if (pixel_opts.fix_factor != 2) {
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

				if (!stop_condition) {
					// must transpose in place, as spectral_matrix.transpose() is const
					spectral_matrix.transposeInPlace();
#ifdef MPI_FOUND
					if (world.size() != 1) {
						solveOneLoop_MPI(
								_data.transpose(),
								pixel_matrix.transpose(),
								spectral_matrix,
								master,
								pixel_opts,
								stop_condition,
								info);
					} else {
#else /* MPI_FOUND */
					{
#endif /* MPI_FOUND */
#ifdef OPENMP_FOUND
						solveOneLoop_OpenMP(
#else /* OPENMP_FOUND */
						solveOneLoop_sequentially(
#endif /* OPENMP_FOUND */
								_data.transpose(),
								pixel_matrix.transpose(),
								spectral_matrix,
								pixel_opts,
								stop_condition,
								auxiliary_constraints,
								info);
					}
					spectral_matrix.transposeInPlace();

					// check criterion
					{
						residual =
								detail::calculateResidual(_data, spectral_matrix, pixel_matrix);
						BOOST_LOG_TRIVIAL(info)
							<< "#" << loop_nr << " 2/2, residual is " << residual;
					}
				}
			}

			// remove ambiguity
			if (spectral_opts.fix_factor == 0) {
				MatrixProductEqualizer equalizer;
				const double scaling_change =
						equalizer(spectral_matrix, pixel_matrix);
				info.replace(IterationInformation::LoopTable, "scaling_change", scaling_change);
				BOOST_LOG_TRIVIAL(info)
					<< "Scaling factor is " << scaling_change;
			}

			if (pixel_opts.fix_factor != 2) {
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
		if (loop_nr > spectral_opts.max_loops)
			BOOST_LOG_TRIVIAL(error)
				<< "Maximum number of loops " << spectral_opts.max_loops
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
#endif /* MPI_FOUND */

		// renormalize both factors
		MatrixProductRenormalizer renormalizer;
		renormalizer(spectral_matrix, pixel_matrix);
	}
}
