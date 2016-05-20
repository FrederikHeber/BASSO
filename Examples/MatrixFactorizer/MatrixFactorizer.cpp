/*
 * MatrixFactorizer.cpp
 *
 *  Created on: Apr 7, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include <Eigen/Dense>
#include <fstream>
#include <string>

#ifdef MPI_FOUND
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;
#endif

#include <boost/chrono.hpp>

#include "MatrixFactorizer/Database/IterationInformation.hpp"
#include "MatrixFactorizer/Helpers/detail.hpp"
#include "MatrixFactorizer/Options/MatrixFactorizerOptions.hpp"
#include "MatrixFactorizer/Solvers/MatrixFactorization.hpp"
#include "MatrixFactorizer/Work/Master.hpp"
#include "MatrixFactorizer/Work/Slave.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "MatrixIO/MatrixIOExceptions.hpp"

int outputSolution(
		const MatrixFactorizerOptions &_opts,
		const Eigen::MatrixXd &_data,
		const Eigen::MatrixXd &_spectral_matrix,
		const Eigen::MatrixXd &_pixel_matrix
		)
{
	int returnstatus = 0;
	/// output solution
	if (!_opts.solution_factor_one_file.string().empty())
		if (!MatrixIO::store(
				_opts.solution_factor_one_file.string(),
				"spectral matrix",
				_spectral_matrix)) {
			returnstatus = 255;
			BOOST_LOG_TRIVIAL(error) <<
					"Failed to write first solution factor file.";
		}
	if (!_opts.solution_factor_two_file.string().empty())
		if (!MatrixIO::store(
				_opts.solution_factor_two_file.string(),
				"pixel matrix",
				_pixel_matrix)) {
			returnstatus = 255;
			BOOST_LOG_TRIVIAL(error) <<
					"Failed to write second solution factor file.";
		}
	const Eigen::MatrixXd product_matrix = _spectral_matrix * _pixel_matrix;
	if (!_opts.solution_product_file.string().empty())
		if (!MatrixIO::store(
				_opts.solution_product_file.string(),
				"solution product",
				product_matrix)) {
			returnstatus = 255;
			BOOST_LOG_TRIVIAL(error) <<
					"Failed to write solution product file.";
 		}
	const Eigen::MatrixXd difference_matrix = _data - product_matrix;
	if (!_opts.solution_difference_file.string().empty())
		if (!MatrixIO::store(
				_opts.solution_difference_file.string(),
				"solution difference",
				difference_matrix)) {
			returnstatus = 255;
			BOOST_LOG_TRIVIAL(error) <<
					"Failed to write solution difference file.";
 		}

	LOG(debug, "Resulting first factor transposed is\n" << _spectral_matrix.transpose());
	LOG(debug, "Resulting second factor is\n" << _pixel_matrix);

	if ((_data.innerSize() <= 10) && (_data.outerSize() <= 10)) {
		LOG(debug, "Data matrix was\n" << _data);
		LOG(debug, "Product matrix is\n" << product_matrix);
		LOG(info, "Difference matrix is\n" << difference_matrix);
	}
	LOG(info, "Norm of difference is " << (_data - product_matrix).norm());

	return returnstatus;
}

int main(int argc, char **argv)
{
	/// start MPI
#ifdef MPI_FOUND
	  mpi::environment env(argc, argv);
	  mpi::communicator world;
#endif /* MPI_FOUND */

#ifdef MPI_FOUND
	if (world.rank() == 0) {
#endif /* MPI_FOUND */
		// show program information
		showVersion(std::string(argv[0]));
		showCopyright();

#ifdef MPI_FOUND
	}
#endif /* MPI_FOUND */

	int returnstatus = 0;
#ifdef MPI_FOUND
	if (world.rank() == 0) {
		LOG(info, "We have one master to distribute and " << (world.size()-1) << " slaves to work on the problem.");
#else /* MPI_FOUND */
#ifdef OPENMP_FOUND
		LOG(info, "Solving with parallel threads.");
#else /* OPENMP_FOUND */
		LOG(info, "A single process solves the problem.");
#endif /* OPENMP_FOUND */
#endif /* MPI_FOUND */
		/// starting timing
		boost::chrono::high_resolution_clock::time_point timing_start =
				boost::chrono::high_resolution_clock::now();

		/// some required parameters
		MatrixFactorizerOptions opts;
		if (returnstatus == 0)
			returnstatus = detail::parseOptions(argc, argv, opts);

		/// parse the matrices

		Eigen::MatrixXd data;
		if (returnstatus == 0) {
			if (opts.sparse) {
				Eigen::SparseMatrix<double, Eigen::RowMajor> sparse_data;
				returnstatus =
						detail::parseDataFile< Eigen::SparseMatrix<double, Eigen::RowMajor> >(
								opts.data_file.string(),
								sparse_data);
				data = sparse_data;
			} else
				returnstatus =
						detail::parseDataFile<Eigen::MatrixXd>(
								opts.data_file.string(),
								data);
		} else {
			LOG(error, "There was an error with the options, exiting.");
		}
		// print parsed matrix and vector if small or high verbosity requested
		if ((data.innerSize() > 10) || (data.outerSize() > 10)) {
			LOG(trace, "We solve for Y=K*X with Y =\n" << data << "." << std::endl);
		} else {
			LOG(info, "We solve for Y=K*X with Y =\n" << data << "." << std::endl);
		}

		if (returnstatus == 0) {

			/// create Database
			IterationInformation info(opts, data.innerSize(), data.outerSize());

			/// perform factorization
			MatrixFactorization factorizer(opts, info
#ifdef MPI_FOUND
				, world
#endif /* MPI_FOUND */
				);
			factorizer(data, returnstatus);

#ifdef MPI_FOUND
			// exchange return status to tell clients whether everything is ok
			// or whether we need to stop execution
			mpi::broadcast(world, returnstatus, 0);
#endif /* MPI_FOUND */

			returnstatus = outputSolution(
					opts,
					data,
					factorizer.spectral_matrix,
					factorizer.pixel_matrix);
		} else {
#ifdef MPI_FOUND
			Master master(world, opts.overall_keys);

			// So far, Slaves are present and expect initial go (or not full_terminate
			// signal). Hence, if something has gone wrong, then we need at least to
			// tell them here before the actual solver loop
			master.sendTerminate();
#else
			LOG(error, "There was an error parsing data matrix, exiting.");
#endif /* MPI_FOUND */
		}

		boost::chrono::high_resolution_clock::time_point timing_end =
				boost::chrono::high_resolution_clock::now();
		LOG(info, "The operation took "
				<< boost::chrono::duration<double>(timing_end - timing_start) << ".");

#ifdef MPI_FOUND
	} else {
		// enter in solve() loop
		Slave worker(world);
		worker();

		// exchange return status
		mpi::broadcast(world, returnstatus, 0);
	}
#endif /* MPI_FOUND */

	if (returnstatus != 0) {
		LOG(error, "There was an error performing the factorization, exiting.");
	}

	/// exit
	return returnstatus;
}
