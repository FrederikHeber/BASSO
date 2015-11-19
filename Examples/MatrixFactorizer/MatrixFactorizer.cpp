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
	if (!_opts.solution_product_file.string().empty())
		if (!MatrixIO::store(
				_opts.solution_product_file.string(),
				"solution product",
				_spectral_matrix * _pixel_matrix)) {
			returnstatus = 255;
			BOOST_LOG_TRIVIAL(error) <<
					"Failed to write solution product file.";
 		}

	BOOST_LOG_TRIVIAL(debug)
		<< "Resulting first factor transposed is\n" << _spectral_matrix.transpose();
	BOOST_LOG_TRIVIAL(debug)
		<< "Resulting second factor is\n" << _pixel_matrix;

	const Eigen::MatrixXd product_matrix = _spectral_matrix * _pixel_matrix;
	if ((_data.innerSize() <= 10) && (_data.outerSize() <= 10)) {
		BOOST_LOG_TRIVIAL(debug)
			<< "Data matrix was\n" << _data;
		BOOST_LOG_TRIVIAL(debug)
			<< "Product matrix is\n" << product_matrix;
		BOOST_LOG_TRIVIAL(info)
			<< "Difference matrix is\n" << _data - product_matrix;
	}
	BOOST_LOG_TRIVIAL(info)
		<< "Norm of difference is " << (_data - product_matrix).norm();

	return returnstatus;
}

int main(int argc, char **argv)
{
	/// start MPI
#ifdef MPI_FOUND
	  mpi::environment env(argc, argv);
	  mpi::communicator world;
#endif

	int returnstatus = 0;
#ifdef MPI_FOUND
	if (world.rank() == 0) {
		BOOST_LOG_TRIVIAL(info)
				<< "We have one master to distribute and "
				<< (world.size()-1) << " slaves to work on the problem.";
#else
		BOOST_LOG_TRIVIAL(info)
			<< "A single process solves the problem.";
#endif
		/// some required parameters
		MatrixFactorizerOptions opts;
		if (returnstatus == 0)
			returnstatus = detail::parseOptions(argc, argv, opts);

		/// parse the matrices
		Eigen::MatrixXd data;
		if (returnstatus == 0)
			returnstatus = detail::parseDataFile(opts.data_file.string(), data);

		/// create Database
		IterationInformation info(opts, data.innerSize(), data.outerSize());

		/// perform factorization
		MatrixFactorization factorizer(opts, info
#ifdef MPI_FOUND
				, world
#endif
				);
		factorizer(data, returnstatus);

#ifdef MPI_FOUND
		// exchange return status to tell clients whether everything is ok
		// or whether we need to stop execution
		mpi::broadcast(world, returnstatus, 0);
#endif

		if (returnstatus == 0)
			returnstatus = outputSolution(
					opts,
					data,
					factorizer.spectral_matrix,
					factorizer.pixel_matrix);

#ifdef MPI_FOUND
	} else {
		// enter in solve() loop
		Slave worker(world);
		worker();

		// exchange return status
		mpi::broadcast(world, returnstatus, 0);
	}
#endif

	/// exit
	return returnstatus;
}
