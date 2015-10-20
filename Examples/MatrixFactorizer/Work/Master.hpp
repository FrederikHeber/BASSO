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

#include "MatrixFactorizer/Work/WorkResult.hpp"

/** This class contains all code related to the master in the master/slave
 * parallelization implementation
 */
class Master
{
public:
	/** Constructor for class Master.
	 *
	 * @param _world world communicator, stored as ref in class.
	 */
	Master(boost::mpi::communicator &_world);

	/** Solve function of Master that essentially just distributes the
	 * data to the Slaves.
	 *
	 *  We go through all columns of \a _rhs or \a _solution, respectively,
	 *  and solve the inverse problem with \a _matrix.
	 *
	 * @param _matrix inverse problem matrix
	 * @param _rhs right-hand side matrix
	 * @param _solution on return solution matrix
	 * @param _loop_nr current outer iteration step (deciding on
	 * 		 when to enforce non-negativity constraint)
	 * @return true - success
	 */
	bool solve(
			const Eigen::MatrixXd &_matrix,
			const Eigen::MatrixXd &_rhs,
			Eigen::MatrixXd &_solution,
			const unsigned int _loop_nr
			);

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
};

#endif /* MPI_FOUND */

#endif /* MATRIXFACTORIZERBASE_WORK_MASTER_HPP_ */
