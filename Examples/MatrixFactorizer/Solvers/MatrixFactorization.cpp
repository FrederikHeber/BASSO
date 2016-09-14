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
	//
	// bullocks: So far, functional analysis has not been very helpful.
	// As SESOP is generally very slow if Y is not an l_2 space, we
	// use the desired spaces for the minimum-norm setting and the l_2
	// space for the data fidelity constraint

	// note on bullocks: because of the factoring the l_2 space would
	// always be the space X, and we could only choose the space Y which
	// is not very helpful,
	// namely: A : X -> Y and  A=KX => X: X -> l_2, K: l_2 -> Y but we
	// calculate X^t, i.e. X^t: l_2 -> X^\ast. Hence, source space is
	// always an l_2
	const_cast<MatrixFactorizerOptions &>(pixel_opts).type_spacey =
			"lp";
	const_cast<MatrixFactorizerOptions &>(pixel_opts).py = 2.;
	const_cast<MatrixFactorizerOptions &>(pixel_opts).powery = 2.;
	const_cast<MatrixFactorizerOptions &>(spectral_opts).type_spacex =
			spectral_opts.type_spacey;
	const_cast<MatrixFactorizerOptions &>(spectral_opts).px = spectral_opts.py;
	const_cast<MatrixFactorizerOptions &>(spectral_opts).powerx = spectral_opts.powery;
	const_cast<MatrixFactorizerOptions &>(spectral_opts).type_spacey =
			"lp";
	const_cast<MatrixFactorizerOptions &>(spectral_opts).py = 2.;
	const_cast<MatrixFactorizerOptions &>(spectral_opts).powery = 2.;
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
	if (!solver_ok) {
		LOG(error, "The minimizer for InRangeSolver could not finish.");
	}
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
	// init solvers and return values
	//!> counter counts the number of successful minimizations, should equal number of columns
	int counter = 0;
	// cannot use std::vector<> here as this requires proper copy and assignment
	// operators which is very complicated due to the internal databases and
	// references for quick and convienient access
	InRangeSolver* solvers[omp_get_max_threads()];
	for (int i=0;i<omp_get_max_threads(); ++i)
		solvers[i] = new InRangeSolver(
				opts,
				opts.overall_keys,
				opts.projection_delta);
	std::vector<AccumulatedValues> projector_values(omp_get_max_threads());
	std::vector<AccumulatedValues> solver_values(omp_get_max_threads());
	if (!stop_condition) {
#pragma omp parallel shared(solvers, projector_values, solver_values) reduction (+ : counter)
		{
			const size_t thread_id = omp_get_thread_num();
			Eigen::VectorXd solution;
#pragma omp for schedule(static)
			for (unsigned int dim = 0; dim < _data.cols(); ++dim) {
				if ((*solvers[thread_id])(fixed_factor,
							_data.col(dim),
							variable_factor.col(dim),
							solution,
							dim,
							auxiliary_constraints
							)) {
					counter = counter+1;
//					#pragma omp critical (solution)
					variable_factor.col(dim) = solution;
				} else {
					LOG(error, "The minimizer for InRangeSolver could not finish on column " << dim << ".");
				}
			}
			projector_values[thread_id] = (*solvers[thread_id]).getAccumulatedProjectorValues();
			solver_values[thread_id] = (*solvers[thread_id]).getAccumulatedProjectorValues();
		} /* end of omp parallel section */

		// gather all values to accumulated in a single solver's databases
		(*solvers[0]).insertProjectorValues(projector_values.begin(), projector_values.end());
		(*solvers[0]).insertSolverValues(solver_values.begin(), solver_values.end());
		// place accumulated values in loop table
		(*solvers[0]).insertAccumulatedProjectorValues(
				info.getLoopTable(), "_projection");
		(*solvers[0]).insertAccumulatedSolverValues(
				info.getLoopTable(), "_minimization");
	}
	stop_condition |= counter != _data.cols();

	// free solvers
	for (int i=0;i<omp_get_max_threads(); ++i)
		delete solvers[i];
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
			LOG(error, "The minimizer for InRangeSolver could not finish on column " << dim << ".");
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
		LOG(info, "#" << loop_nr << " starting residual is " << residual);
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
				LOG(error, "Parsing of either or both factors failed, aborting.");
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
			LOG(info, "Setting spectral matrix K to random starting values.");
			spectral_matrix = Eigen::MatrixXd(_data.rows(), spectral_opts.sparse_dim);
			detail::constructRandomMatrix(spectral_matrix);
		} else {
			if ((spectral_matrix.rows() == _data.rows())
					&& (spectral_matrix.cols() == spectral_opts.sparse_dim)) {
				LOG(info, "Using spectral matrix K parsed from file.");
			} else {
				LOG(error, "Parsed spectral matrix has dimensions "
						<< spectral_matrix.rows() << "," << spectral_matrix.cols() << ", while expecting " << _data.rows() << "," << spectral_opts.sparse_dim);
				_returnstatus = 255;
#ifdef MPI_FOUND
				master.sendTerminate();
#endif /* MPI_FOUND */
				return;
			}
		}
		if (!parse_status) {
			LOG(info, "Setting pixel matrix X to zero.");
			pixel_matrix = Eigen::MatrixXd(spectral_opts.sparse_dim, _data.cols());
			detail::constructZeroMatrix(
					pixel_matrix);
		} else {
			if ((pixel_matrix.rows() == spectral_opts.sparse_dim)
					&& (pixel_matrix.cols() == _data.cols())) {
				LOG(info, "Using pixel matrix X parsed from file.");
			} else {
				LOG(error, "Parsed pixel matrix has dimensions "
						<< pixel_matrix.rows() << "," << pixel_matrix.cols() << ", while expecting " << spectral_opts.sparse_dim << "," <<  _data.cols());
				_returnstatus = 255;
#ifdef MPI_FOUND
				master.sendTerminate();
#endif /* MPI_FOUND */
				return;
			}
		}

		// construct transposed version of the data matrix (ColMajor is default)
		const Eigen::MatrixXd data_transposed = _data.transpose();

		// create outer ("loop") stopping criteria
		StoppingArguments stopping_args;
		stopping_args.setTolerance(spectral_opts.residual_threshold);
		stopping_args.setMaxIterations(spectral_opts.max_loops);
		StoppingCriteriaFactory stop_factory;
		StoppingCriterion::ptr_t stopping_criterion =
				stop_factory.create(spectral_opts.factorization_stopping_criteria, stopping_args);

		// create auxiliary constraints
		AuxiliaryConstraintsFactory constraint_factory;
		AuxiliaryConstraints::ptr_t auxiliary_constraints =
				constraint_factory.create(spectral_opts.auxiliary_constraints);

		/// create feasible solution starting points
		if (auxiliary_constraints) {
			double norm_difference = spectral_matrix.norm();
			if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
				LOG(trace, "Spectral matrix K^t before constraining is \n" << spectral_matrix.transpose());
			} else {
				LOG(info, "Spectral matrix K^t before constraining is \n" << spectral_matrix.transpose());
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
//				LOG(info, "tempvector of col " << i << " before constraining is " << *tempvector);
				(*auxiliary_constraints)(tempvector);
//				LOG(info, "tempvector of col " << i << " after constraining is " << *tempvector);
				const Eigen::Ref<Eigen::MatrixXd::ColXpr> col = spectral_matrix.col(i);
				VectorSetter::set(col, tempvector);
			}
			spectral_matrix.transposeInPlace();
			norm_difference -= spectral_matrix.norm();
			LOG(info, "Spectral matrix K^t's l2 norm due to constraints changed by " << norm_difference);

			if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
				LOG(trace, "Spectral matrix K^t after constraining is \n" << spectral_matrix.transpose());
			} else {
				LOG(info, "Spectral matrix K^t after constraining is \n" << spectral_matrix.transpose());
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
				LOG(debug, "======================== #" << loop_nr << "/1 ==================");

		//		renormalizeMatrixByTrace(spectral_matrix);

				if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
					LOG(trace, "Current spectral matrix K^t is\n" << spectral_matrix.transpose());
				} else {
					LOG(info, "Current spectral matrix K^t is\n" << spectral_matrix.transpose());
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
					LOG(trace, "Resulting pixel matrix X is\n" << pixel_matrix);
				} else {
					LOG(info, "Resulting pixel matrix X is\n" << pixel_matrix);
				}

				// check criterion
				{
					residual =
							detail::calculateResidual(_data, spectral_matrix, pixel_matrix);
					info.replace(IterationInformation::LoopTable, "residual", residual);
					LOG(info, "#" << loop_nr << " 1/2, residual is " << residual);
				}

				// submit loop tuple
				info.addTuple(IterationInformation::LoopTable);
			}

			if (pixel_opts.fix_factor != 2) {
				LOG(debug, "======================== #" << loop_nr << "/2 ==================");

				info.replace(IterationInformation::LoopTable, "loop_nr", 2*(int)(loop_nr));

		//		renormalizeMatrixByTrace(pixel_matrix);

				if ((pixel_matrix.innerSize() > 10) || (pixel_matrix.outerSize() > 10)) {
					LOG(trace, "Current pixel matrix X is\n" << pixel_matrix);
				} else {
					LOG(info, "Current pixel matrix X is\n" << pixel_matrix);
				}

				if (!stop_condition) {
					// must transpose in place, as spectral_matrix.transpose() is const
					spectral_matrix.transposeInPlace();
#ifdef MPI_FOUND
					if (world.size() != 1) {
						solveOneLoop_MPI(
								data_transposed,
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
								data_transposed,
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
						LOG(info, "#" << loop_nr << " 2/2, residual is " << residual);
					}
				}
			}

			// remove ambiguity
			if (spectral_opts.fix_factor == 0) {
				MatrixProductEqualizer equalizer;
				const double scaling_change =
						equalizer(spectral_matrix, pixel_matrix);
				info.replace(IterationInformation::LoopTable, "scaling_change", scaling_change);
				LOG(info, "Scaling factor is " << scaling_change);
			}

			if (pixel_opts.fix_factor != 2) {
				if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
					LOG(trace, "Resulting spectral K^t matrix is\n" << spectral_matrix.transpose());
				} else {
					LOG(info, "Resulting spectral K^t matrix is\n" << spectral_matrix.transpose());
				}

				// check criterion
				{
					residual =
							detail::calculateResidual(_data, spectral_matrix, pixel_matrix);
					LOG(info, "#" << loop_nr << ", residual is " << residual);
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
		LOG(info, "Iteration was stopped by "
				<< stopping_criterion->whoIsTrue(
						boost::chrono::duration<double>(0),
						loop_nr,
						residual,
						1.));

		if (loop_nr > spectral_opts.max_loops) {
			LOG(error, "Maximum number of loops " << spectral_opts.max_loops << " exceeded, stopping iteration.");
		} else {
			LOG(info, "Loop iteration performed " << loop_nr << " times.");
		}

		/// end timing
		boost::chrono::high_resolution_clock::time_point timing_end =
				boost::chrono::high_resolution_clock::now();
		LOG(info, "The operation took "
				<< boost::chrono::duration<double>(timing_end - timing_start).count() << " seconds.");

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
