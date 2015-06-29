/*
 * MatrixFactorization.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_SOLVERS_MATRIXFACTORIZATION_HPP_
#define MATRIXFACTORIZERBASE_SOLVERS_MATRIXFACTORIZATION_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#ifdef MPI_FOUND
#include <boost/mpi/communicator.hpp>
#endif

class MatrixFactorizerOptions;
class IterationInformation;

/** Functor that performs a Matrix Factorization.
 *
 */
struct MatrixFactorization
{
	/** Constructor for class MatrixFactorization.
	 *
	 * @param opts options for controlling inverse problem solving
	 * @param info database connection
	 */
	MatrixFactorization(
			const MatrixFactorizerOptions &_opts,
			IterationInformation &_info
#ifdef MPI_FOUND
			, boost::mpi::communicator &_world
#endif
			);

	/** Performs the factorization.
	 *
	 * @param _data matrix to factorizer
	 * @param _returnstatus indicates whether to stop right away or
	 *        attempt solving, on end contains return status: 0 - success
	 */
	void operator()(
			const Eigen::MatrixXd &_data,
			int &_returnstatus
			);

private:
	//!> options for solving inverse problems
	const MatrixFactorizerOptions &opts;
	//!> database connection
	IterationInformation &info;
#ifdef MPI_FOUND
		boost::mpi::communicator &world;
#endif
};



#endif /* MATRIXFACTORIZERBASE_SOLVERS_MATRIXFACTORIZATION_HPP_ */
