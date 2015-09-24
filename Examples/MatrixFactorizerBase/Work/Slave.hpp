/*
 * Slave.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_WORK_SLAVE_HPP_
#define MATRIXFACTORIZERBASE_WORK_SLAVE_HPP_

#include "BassoConfig.h"

#ifdef MPI_FOUND
#include <boost/mpi/communicator.hpp>

/** This class implements the Slave that solves specific inverse problems
 * sent to him by the Master.
 */
class Slave
{
public:
	/** Constructor for class Master.
	 *
	 * @param _world world communicator, stored as ref in class.
	 */
	Slave(boost::mpi::communicator &_world);

	/** Worker function for arbitrary "slave" process to work on inverse
	 * problems for solving matrix factorization problem.
	 *
	 * Workers operate in two nested loops:
	 * -# the inner loop is for handling the columns of a single expression
	 *    \f$ Y_i = K*X_i \f$ or $\f$ Y^t_i = X^t * K^t_i \f$, respectively.
	 * -# the outer loop is for the alternating between the two matrix factors
	 *    K and X. It resets the inner loop with respect to the data that is
	 *    "global" therein.
	 *
	 */
	void operator()();

private:
	//!> reference to world communicator
	boost::mpi::communicator &world;
};

#endif /* MPI_FOUND */

#endif /* MATRIXFACTORIZERBASE_WORK_SLAVE_HPP_ */
