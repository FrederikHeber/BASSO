/*
 * InverseProblemSolver.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef SOLVERS_INVERSEPROBLEMSOLVER_HPP_
#define SOLVERS_INVERSEPROBLEMSOLVER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/types.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"

class CommandLineOptions;

/** Functor for solving inverse problems
 *
 */
struct InverseProblemSolver
{
	/** Constructor for class InverseProblemSolver.
	 *
	 * @param _database ref to database
	 * @param _opts options defining manner of inverse problem solving
	 * @param _checkTrueSolution whether true solution by SVD should be
	 * 			calculated, only possible for l2 norms
	 */
	InverseProblemSolver(
			Database_ptr_t &_database,
			const CommandLineOptions &_opts,
			const bool _checkTrueSolution
			);

	/** Functor that solves the inverse problem.
	 *
	 * @param _inverseproblem structure defining the inverse problem
	 * @param _startingvalue starting value
	 * @return true - success
	 */
	GeneralMinimizer::ReturnValues operator()(
			InverseProblem_ptr_t &_inverseproblem,
			const Eigen::VectorXd &_startingvalue
			);

	/** Functor that solves the inverse problem.
	 *
	 * @param _matrix problem matrix
	 * @param _rhs right-hand side
	 * @param _startingvalue starting value
	 * @param _solution to contain solution on exit
	 * @return true - success, false - else
	 */
	bool operator()(
			const Eigen::MatrixXd &_matrix,
			const Eigen::MatrixXd &_rhs,
			const Eigen::VectorXd &_startingvalue,
			Eigen::VectorXd &_solution
			);
private:
	//!> database reference for storing iteration information
	Database_ptr_t database;
	//!> storing options on the manner of solving
	const CommandLineOptions &opts;
	//!> whether to calculate true solution in l2 case from SVD
	const bool checkTrueSolution;
};



#endif /* SOLVERS_INVERSEPROBLEMSOLVER_HPP_ */
