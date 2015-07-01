/*
 * Slave.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Slave.hpp"

#ifdef MPI_FOUND
#include <boost/mpi/environment.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/nonblocking.hpp>

namespace mpi = boost::mpi;

#include "Log/Logging.hpp"
#include "MatrixFactorizerBase/Helpers/detail.hpp"
#include "MatrixFactorizerBase/Options/MatrixFactorizerOptions.hpp"
#include "MatrixFactorizerBase/Solvers/InRangeSolver.hpp"
#include "MatrixFactorizerBase/Work/WorkPackage.hpp"
#include "MatrixFactorizerBase/Work/WorkResult.hpp"
#include "Minimizations/Elements/Eigen_matrix_serialization.hpp"

Slave::Slave(boost::mpi::communicator &_world) :
		world(_world)
{}

void Slave::operator()()
{
	// get global information
	BOOST_LOG_TRIVIAL(debug)
			<< "#" << world.rank() << " - getting options.";
	MatrixFactorizerOptions opts;
	mpi::broadcast(world, opts, 0);
	InRangeSolver solver(opts);

	bool full_terminate = false;
	while (!full_terminate) {
		// check whether we have to terminate
		mpi::broadcast(world, full_terminate, 0);
		if (full_terminate)
			break;

		// send id
		BOOST_LOG_TRIVIAL(debug)
				<< "#" << world.rank() << " - initiating by sending id.";
		world.send(0, detail::InitialId, world.rank());

		// get global information
		BOOST_LOG_TRIVIAL(debug)
				<< "#" << world.rank() << " - getting matrix.";
		Eigen::MatrixXd matrix;
		mpi::broadcast(world, matrix, 0);
		unsigned int loop_nr = 0;
		mpi::broadcast(world, loop_nr, 0);

		if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
			BOOST_LOG_TRIVIAL(trace)
					<< "#" << world.rank() << " - got global data\n"
					<< matrix;
		} else {
			BOOST_LOG_TRIVIAL(debug)
					<< "#" << world.rank() << " - got global data\n"
					<< matrix;
		}

		// we stop working only when we get the termination signal
		// from the master
		bool terminate = false;
		while ((!terminate) && (!full_terminate)) {
			/// wait for receiving data
			WorkPackage package;
			BOOST_LOG_TRIVIAL(debug)
					<< "#" << world.rank() << " - receiving next work package.";
			world.recv(0, detail::ColumnWise, package);
			const int col = package.col;
			terminate |= (col == -1);
			if (!terminate) {
				const Eigen::VectorXd &rhs = package.rhs;
				const Eigen::VectorXd &solution_startvalue =
						package.solution_startvalue;

				if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
					BOOST_LOG_TRIVIAL(trace)
							<< "#" << world.rank() << " got problem rhs, col "
							<< col << "\n" << rhs.transpose();
				} else {
					BOOST_LOG_TRIVIAL(debug)
							<< "#" << world.rank() << " got problem rhs, col "
							<< col << "\n" << rhs.transpose();
				}

				/// work on data
				BOOST_LOG_TRIVIAL(debug)
						<< "#" << world.rank() << " - working.";
				Eigen::VectorXd solution;
				const bool solve_ok =
						solver(matrix,
								rhs,
								solution_startvalue,
								solution,
								col,
								loop_nr);

				if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
					BOOST_LOG_TRIVIAL(trace)
							<< "#" << world.rank() << " sending solution, col "
							<< col << "\n" << solution.transpose();
				} else {
					BOOST_LOG_TRIVIAL(debug)
							<< "#" << world.rank() << " sending solution, col "
							<< col << "\n" << solution.transpose();
				}

				/// return result
				BOOST_LOG_TRIVIAL(debug)
						<< "#" << world.rank() << " - sending solution.";
				WorkResult result(solution, solve_ok, col);
				world.send(0, detail::ColumnWise, result);
			} else {
				BOOST_LOG_TRIVIAL(debug)
						<< "#" << world.rank() << " - terminating.";
			}
		}
	}
}

#endif /* MPI_FOUND */