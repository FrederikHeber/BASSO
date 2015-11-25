/*
 * Master.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_WORK_MASTER_HPP_
#define MATRIXFACTORIZERBASE_WORK_MASTER_HPP_

#include "BassoConfig.h"

#ifdef MPI_FOUND
#include <boost/mpi/communicator.hpp>

#include <vector>

#include <Eigen/Dense>

#include "MatrixFactorizer/Database/InnerProblemDatabase.hpp"
#include "MatrixFactorizer/Work/WorkResult.hpp"

class CommandLineOptions;
class Table;

/** This class contains all code related to the master in the master/slave
 * parallelization implementation
 */
class Master
{
public:
	/** Constructor for class Master.
	 *
	 * @param _world world communicator, stored as ref in class.
	 * @param _overall_keys set of keys to accumulate from solver's overall
	 * 			table
	 */
	Master(
			boost::mpi::communicator &_world,
			const InnerProblemDatabase::keys_t &_overall_keys);

	/** Solve function of Master that essentially just distributes the
	 * data to the Slaves.
	 *
	 *  We go through all columns of \a _rhs or \a _solution, respectively,
	 *  and solve the inverse problem with \a _matrix.
	 *
	 * @param _opts options to steer how inverse problem is solved
	 * @param _matrix inverse problem matrix
	 * @param _rhs right-hand side matrix
	 * @param _solution on return solution matrix
	 * @return true - success
	 */
	bool solve(
			const CommandLineOptions &_opts,
			const Eigen::MatrixXd &_matrix,
			const Eigen::MatrixXd &_rhs,
			Eigen::MatrixXd &_solution
			);

	/** Insert accumulated and gathered values from projector problem into
	 * given \a _table.
	 *
	 * @param _table table to insert values to
	 */
	void insertAccumulatedProjectorValues(
			Table &_table) const;

	/** Insert accumulated and gathered values from solver problem into
	 * given \a _table.
	 *
	 * @param _table table to insert values to
	 */
	void insertAccumulatedSolverValues(
			Table &_table) const;

	/** Resets the internal databases for projector problem
	 *
	 */
	void resetAccumulatedProjectorValues();

	/** Resets the internal databases for solver problem.
	 *
	 */
	void resetAccumulatedSolverValues();

	/** Sends the terminate signal to the Slaves.
	 *
	 */
	void sendTerminate();

private:
	/** Handles results obtained asynchronously.
	 *
	 * \param _result asynchronous result
	 * \param _results results array where the \a result's contents is stored
	 * \param _solution solution vector where obtained solution is stored
	 */
	bool handleResult(
			const boost::mpi::status &_result,
			const std::vector<WorkResult> &_results,
			Eigen::MatrixXd &_solution
			);
private:
	//!> reference to world communicator
	boost::mpi::communicator &world;
	//!> ref to overall keys to broadcast to clients
	const InnerProblemDatabase::keys_t &overall_keys;

	//!> table database for projector iteration information
	Database::Database_ptr_t projector_db;
	//!> table database for solver iteration information
	Database::Database_ptr_t solver_db;
};

#endif /* MPI_FOUND */

#endif /* MATRIXFACTORIZERBASE_WORK_MASTER_HPP_ */
