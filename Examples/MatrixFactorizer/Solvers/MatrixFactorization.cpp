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
#include "MatrixFactorizer/IterationChecks/RelativeResidualChecker.hpp"
#include "MatrixFactorizer/IterationChecks/ResidualChecker.hpp"
#include "MatrixFactorizer/IterationChecks/MaxIterationsChecker.hpp"
#include "MatrixFactorizer/ScalingAmbiguityRemover/MatrixProductEqualizer.hpp"
#include "MatrixFactorizer/ScalingAmbiguityRemover/MatrixProductRenormalizer.hpp"
#include "MatrixFactorizer/Solvers/InRangeSolver.hpp"

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

	if (_returnstatus != 0) {
		master.sendTerminate();
	} else
#endif
	{
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
		MaxIterationsChecker MaxIterations(opts.max_loops);
		RelativeChangeResidualChecker RelativeResidual(opts.residual_threshold);
		ResidualChecker Residual(opts.residual_threshold);

		bool stop_condition =
				RelativeResidual(residual)
				|| Residual(residual)
				|| MaxIterations(loop_nr);

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
									dim,
									loop_nr);
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
								pixel_matrix,
								loop_nr);
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
									dim,
									loop_nr);
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
								spectral_matrix,
								loop_nr);
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
				stop_condition =
						RelativeResidual(residual)
						|| Residual(residual)
						|| MaxIterations(loop_nr);
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
		info.replace(IterationInformation::OverallTable, "loops", (int)loop_nr);
		info.replace(IterationInformation::OverallTable, "residual", residual);
		info.addTuple(IterationInformation::OverallTable);

#ifdef MPI_FOUND
		master.sendTerminate();
#endif

		// renormalize both factors
		MatrixProductRenormalizer renormalizer;
		renormalizer(spectral_matrix, pixel_matrix);

	}
}
