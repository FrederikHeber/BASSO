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
#endif
		/// starting timing
		boost::chrono::high_resolution_clock::time_point timing_start =
				boost::chrono::high_resolution_clock::now();

		/// some required parameters
		MatrixFactorizerOptions opts;
		if (returnstatus == 0)
			returnstatus = detail::parseOptions(argc, argv, opts);

#ifdef MPI_FOUND
		// send round options
		BOOST_LOG_TRIVIAL(debug)
				<< "#0 - broadcasting options.";
		mpi::broadcast(world, opts, 0);
#endif

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
		// exchange return status
		mpi::broadcast(world, returnstatus, 0);
#endif

		/// finish timing
		boost::chrono::high_resolution_clock::time_point timing_end =
				boost::chrono::high_resolution_clock::now();
		BOOST_LOG_TRIVIAL(info) << "The operation took "
				<< boost::chrono::duration<double>(timing_end - timing_start)
				<< ".";
#ifdef MPI_FOUND
	} else {
		// enter in solve() loop
		Slave worker(world);
		worker();

		// exchange return status
		mpi::broadcast(world, returnstatus, 0);
	}

	// End MPI
	MPI_Finalize ();
#endif

	/// exit
	return returnstatus;
}
