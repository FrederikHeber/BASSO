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
#include "MatrixFactorizer/Helpers/detail.hpp"
#include "MatrixFactorizer/Options/MatrixFactorizerOptions.hpp"
#include "MatrixFactorizer/Solvers/InRangeSolver.hpp"
#include "MatrixFactorizer/Work/WorkPackage.hpp"
#include "MatrixFactorizer/Work/WorkResult.hpp"
#include "Minimizations/Elements/Eigen_matrix_serialization.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraintsFactory.hpp"

Slave::Slave(boost::mpi::communicator &_world) :
		world(_world)
{}

void Slave::operator()()
{
	bool full_terminate = false;
	while (!full_terminate) {
		// check whether we have to terminate
		mpi::broadcast(world, full_terminate, 0);
		if (full_terminate)
			break;

		// send id
		LOG(info, "#" << world.rank() << " - initiating by sending id.");
		world.send(0, detail::InitialId, world.rank());

		// get global information
		LOG(info, "#" << world.rank() << " - getting options.");
		MatrixFactorizerOptions opts;
		mpi::broadcast(world, opts, 0);

		LOG(debug, "#" << world.rank() << " - getting constraints.");
		std::string auxiliary_constraints_string;
		mpi::broadcast(world, auxiliary_constraints_string, 0);

		LOG(debug, "#" << world.rank() << " - getting matrix.");
		Eigen::MatrixXd matrix;
		mpi::broadcast(world, matrix, 0);

		if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
			LOG(trace, "#" << world.rank() << " - got global data\n" << matrix);
		} else {
			LOG(debug, "#" << world.rank() << " - got global data\n" << matrix);
		}

		InRangeSolver solver(opts, opts.overall_keys, opts.projection_delta);

		// create auxiliary constraints
		AuxiliaryConstraintsFactory constraint_factory;
		AuxiliaryConstraints::ptr_t auxiliary_constraints =
				constraint_factory.create(auxiliary_constraints_string);

		// we stop working only when we get the termination signal
		// from the master
		bool terminate = false;
		while ((!terminate) && (!full_terminate)) {
			/// wait for receiving data
			WorkPackage package;
			LOG(debug, "#" << world.rank() << " - receiving next work package.");
			world.recv(0, detail::ColumnWise, package);
			const int col = package.col;
			terminate |= (col == -1);
			if (!terminate) {
				const Eigen::VectorXd &rhs = package.rhs;
				const Eigen::VectorXd &solution_startvalue =
						package.solution_startvalue;

				if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
					LOG(trace, "#" << world.rank() << " got problem rhs, col " << col << "\n" << rhs.transpose());
					LOG(trace, "#" << world.rank() << " got startvalue, col " << col << "\n" << solution_startvalue.transpose());
				} else {
					LOG(debug, "#" << world.rank() << " got problem rhs, col " << col << "\n" << rhs.transpose());
					LOG(debug, "#" << world.rank() << " got startvalue, col " << col << "\n" << solution_startvalue.transpose());
				}

				/// work on data
				LOG(debug, "#" << world.rank() << " - working.");
				Eigen::VectorXd solution;
				const bool solve_ok =
						solver(matrix,
								rhs,
								solution_startvalue,
								solution,
								col,
								auxiliary_constraints);

				if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
					LOG(trace, "#" << world.rank() << " sending solution, col " << col << "\n" << solution.transpose());
				} else {
					LOG(debug, "#" << world.rank() << " sending solution, col " << col << "\n" << solution.transpose());
				}

				/// return result
				LOG(debug, "#" << world.rank() << " - sending solution.");
				WorkResult result(solution, solve_ok, col);
				world.send(0, detail::ColumnWise, result);
			} else {
				LOG(debug, "#" << world.rank() << " - terminating.");
			}
		}
		// gather accumulated values from solved inner problems
		{
			const AccumulatedValues &in_values =
					solver.getAccumulatedProjectorValues();
			LOG(debug, "#" << world.rank() << " sending " << in_values.getNumberOfValues() << " accumulated values from projector problem.");
			boost::mpi::gather(world, in_values, 0);
		}
		{
			const AccumulatedValues &in_values =
					solver.getAccumulatedSolverValues();
			LOG(debug, "#" << world.rank() << " sending " << in_values.getNumberOfValues() << " accumulated values from minimization problem.");
			boost::mpi::gather(world, in_values, 0);
		}
	}
}

#endif /* MPI_FOUND */
