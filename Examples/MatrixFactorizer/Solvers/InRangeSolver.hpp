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

#include "MatrixFactorizer/Database/InnerProblemDatabase.hpp"
#include "Minimizations/types.hpp"

class CommandLineOptions;
class InnerProblemDatabase;

/** Solver combines projection onto range and inverse problem solving with
 * SESOP that requires right-hand side to be in range.
 *
 */
struct InRangeSolver
{
	/** Constructor for class Solver.
	 *
	 * @param _opts options defining the manner of solving.
	 * @param _overall_keys set of keys to accumulate from solver's overall
	 * 			table
	 */
	InRangeSolver(
			const CommandLineOptions &_opts,
			const InnerProblemDatabase::keys_t &_overall_keys
			);

	/** Functor that projects onto range and solves.
	 *
	 * @param _matrix problem matrix
	 * @param _rhs right-hand side
	 * @param _solution_start initial value for solution
	 * @param _solution on return containing solution
	 * @param _dim unique index (for output statements only)
	 * @return true - success
	 */
	bool operator()(
			const Eigen::MatrixXd &_matrix,
			const Eigen::VectorXd &_rhs,
			const Eigen::VectorXd &_solution_start,
			Eigen::VectorXd &_solution,
			const unsigned int _dim
			);

	/** Insert accumulated values from projector problem into given \a _table.
	 *
	 * @param _table table to insert values to
	 */
	void insertAccumulatedProjectorValues(
			Table &_table) const;

	/** Insert accumulated values from solver problem into given \a _table.
	 *
	 * @param _table table to insert values to
	 */
	void insertAccumulatedSolverValues(
			Table &_table) const;

private:
	const CommandLineOptions &opts;
	//!> table database for projector iteration information
	Database_ptr_t projector_db;
	//!> table database for solver iteration information
	Database_ptr_t solver_db;

	InnerProblemDatabase &projectorDB;
	InnerProblemDatabase &solverDB;
};



#endif /* MATRIXFACTORIZERBASE_SOLVERS_INRANGESOLVER_HPP_ */
