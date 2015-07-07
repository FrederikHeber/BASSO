/*
 * Solver.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_SOLVERS_INRANGESOLVER_HPP_
#define MATRIXFACTORIZERBASE_SOLVERS_INRANGESOLVER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/types.hpp"

class CommandLineOptions;

/** Solver combines projection onto range and inverse problem solving with
 * SESOP that requires right-hand side to be in range.
 *
 */
struct InRangeSolver
{
	/** Constructor for class Solver.
	 *
	 * @param _opts options defining the manner of solving.
	 */
	InRangeSolver(const CommandLineOptions &_opts);

	/** Constructor for class Solver.
	 *
	 * @param _opts options defining the manner of solving.
	 * @param _db database to store iteration information
	 */
	InRangeSolver(
			const CommandLineOptions &_opts,
			Database_ptr_t &_db);

	/** Functor that projects onto range and solves.
	 *
	 * @param _matrix problem matrix
	 * @param _rhs right-hand side
	 * @param _solution_start initial value for solution
	 * @param _solution on return containing solution
	 * @param _dim unique index (for output statements only)
	 * @param _loop_nr current outer iteration step (deciding on
	 * 		 when to enforce non-negativity constraint)
	 * @return true - success
	 */
	bool operator()(
			const Eigen::MatrixXd &_matrix,
			const Eigen::VectorXd &_rhs,
			const Eigen::VectorXd &_solution_start,
			Eigen::VectorXd &_solution,
			const unsigned int _dim,
			const unsigned int _loop_nr
			);

private:

	const CommandLineOptions &opts;
	//!> mock database as projector and solver require them
	Database_ptr_t mock_db;
};



#endif /* MATRIXFACTORIZERBASE_SOLVERS_INRANGESOLVER_HPP_ */
