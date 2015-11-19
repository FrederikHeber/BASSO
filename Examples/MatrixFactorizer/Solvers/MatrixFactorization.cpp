/*
 * MatrixFactorization.cpp
 *
 *  Created on: Jun 29, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "MatrixFactorization.hpp"

#include <fstream>
#include <string>

#include "Log/Logging.hpp"
#include "MatrixFactorizer/Database/IterationInformation.hpp"
#include "MatrixFactorizer/Helpers/detail.hpp"
#include "MatrixFactorizer/Options/MatrixFactorizerOptions.hpp"
#include "MatrixFactorizer/Work/Master.hpp"
#include "MatrixFactorizer/ScalingAmbiguityRemover/MatrixProductEqualizer.hpp"
#include "MatrixFactorizer/ScalingAmbiguityRemover/MatrixProductRenormalizer.hpp"
#include "MatrixFactorizer/Solvers/InRangeSolver.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriteriaFactory.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingArguments.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion.hpp"

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
	Master master(world);

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

		/// construct solution starting points
		spectral_matrix = Eigen::MatrixXd(_data.rows(), opts.sparse_dim);
		pixel_matrix = Eigen::MatrixXd(opts.sparse_dim, _data.cols());
		detail::constructStartingMatrices(
				spectral_matrix,
				pixel_matrix,
				opts.type_spacex.find("nonnegative") != std::string::npos);

		/// iterate over the two factors
		double residual =
				detail::calculateResidual(_data, spectral_matrix, pixel_matrix);
		BOOST_LOG_TRIVIAL(info)
			<< "#" << loop_nr << " 1/2, residual is " << residual;
		info.replace(IterationInformation::LoopTable, "residual", residual);

		// create outer ("loop") stopping criteria
		StoppingArguments stopping_args;
		stopping_args.setTolerance(opts.residual_threshold);
		stopping_args.setMaxIterations(opts.max_loops);
		std::string stopping_criteria(
				"MaxIterationCount || RelativeChangeResiduum || Residuum");
		StoppingCriteriaFactory stop_factory;
		StoppingCriterion::ptr_t stopping_criterion =
				stop_factory.create(stopping_criteria, stopping_args);

		bool stop_condition = (*stopping_criterion)(
				boost::chrono::duration<double>(0),
				loop_nr,
				residual,
				1.);

		// submit loop tuple
		info.addTuple(IterationInformation::LoopTable);

		while (!stop_condition) {
			// update loop count
			++loop_nr;
			info.replace(IterationInformation::LoopTable, "loop_nr", (int)loop_nr);

			BOOST_LOG_TRIVIAL(debug)
				<< "======================== #" << loop_nr << "/1 ==================";

	//		renormalizeMatrixByTrace(spectral_matrix);

			if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
				BOOST_LOG_TRIVIAL(trace)
						<< "Current spectral matrix is\n" << spectral_matrix;
			} else {
				BOOST_LOG_TRIVIAL(info)
						<< "Current spectral matrix is\n" << spectral_matrix;
			}

#ifdef MPI_FOUND
			if (world.size() == 1)
#endif
			{
				InRangeSolver solver(spectral_opts);
				for (unsigned int dim = 0; dim < _data.cols(); ++dim) {
					Eigen::VectorXd solution;
					stop_condition &=
							!solver(spectral_matrix,
									_data.col(dim),
									pixel_matrix.col(dim),
									solution,
									dim
									);
					pixel_matrix.col(dim) = solution;
				}
			}
#ifdef MPI_FOUND
			else {
				stop_condition &=
						!master.solve(
								spectral_opts,
								spectral_matrix,
								_data,
								pixel_matrix
								);
			}
#endif

			if ((pixel_matrix.innerSize() > 10) || (pixel_matrix.outerSize() > 10)) {
				BOOST_LOG_TRIVIAL(trace)
						<< "Resulting pixel matrix is\n" << pixel_matrix;
			} else {
				BOOST_LOG_TRIVIAL(info)
						<< "Resulting pixel matrix is\n" << pixel_matrix;
			}

			// check criterion
			{
				residual =
						detail::calculateResidual(_data, spectral_matrix, pixel_matrix);
				BOOST_LOG_TRIVIAL(info)
					<< "#" << loop_nr << " 1/2, residual is " << residual;
			}

			BOOST_LOG_TRIVIAL(debug)
				<< "======================== #" << loop_nr << "/2 ==================";

	//		renormalizeMatrixByTrace(pixel_matrix);

			if ((pixel_matrix.innerSize() > 10) || (pixel_matrix.outerSize() > 10)) {
				BOOST_LOG_TRIVIAL(trace)
						<< "Current pixel matrix is\n" << pixel_matrix;
			} else {
				BOOST_LOG_TRIVIAL(info)
						<< "Current pixel matrix is\n" << pixel_matrix;
			}

			// must transpose in place, as spectral_matrix.transpose() is const
			spectral_matrix.transposeInPlace();
#ifdef MPI_FOUND
			if (world.size() == 1)
# endif
			{
				InRangeSolver solver(pixel_opts);
				for (unsigned int dim = 0; dim < _data.rows(); ++dim) {
					Eigen::VectorXd solution;
					stop_condition &=
							solver(pixel_matrix.transpose(),
									_data.row(dim).transpose(),
									spectral_matrix.col(dim),
									solution,
									dim
									);
					spectral_matrix.col(dim) = solution;
				}
			}
#ifdef MPI_FOUND
			else {
				stop_condition &=
						!master.solve(
								pixel_opts,
								pixel_matrix.transpose(),
								_data.transpose(),
								spectral_matrix
								);
			}
#endif
			spectral_matrix.transposeInPlace();

			// remove ambiguity
			MatrixProductEqualizer equalizer;
			equalizer(spectral_matrix, pixel_matrix);

			if ((spectral_matrix.innerSize() > 10) || (spectral_matrix.outerSize() > 10)) {
				BOOST_LOG_TRIVIAL(trace)
						<< "Resulting spectral matrix is\n" << spectral_matrix;
			} else {
				BOOST_LOG_TRIVIAL(info)
						<< "Resulting spectral matrix is\n" << spectral_matrix;
			}

			// check criterion
			{
				residual =
						detail::calculateResidual(_data, spectral_matrix, pixel_matrix);
				BOOST_LOG_TRIVIAL(info)
					<< "#" << loop_nr << " 2/2, residual is " << residual;
				info.replace(IterationInformation::LoopTable, "residual", residual);
				stop_condition = (*stopping_criterion)(
						boost::chrono::duration<double>(0),
						loop_nr,
						residual,
						1.);
			}

			// submit loop tuple
			info.addTuple(IterationInformation::LoopTable);
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
