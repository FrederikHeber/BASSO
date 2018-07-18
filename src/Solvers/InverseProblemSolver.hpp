/*
 * InverseProblemSolver.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef SOLVERS_INVERSEPROBLEMSOLVER_HPP_
#define SOLVERS_INVERSEPROBLEMSOLVER_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Solvers/FeasibilityProblem.hpp"

class CommandLineOptions;

/** Functor for solving inverse problems
 *
 */
struct InverseProblemSolver : public FeasibilityProblem
{
	/** Constructor for class InverseProblemSolver.
	 *
	 * @param _inverseproblem structure defining the inverse problem
	 * @param _database ref to database
	 * @param _opts options defining manner of inverse problem solving
	 * @param _checkTrueSolution whether true solution by SVD should be
	 * 			calculated, only possible for l2 norms
	 */
	InverseProblemSolver(
			InverseProblem_ptr_t &_inverseproblem,
			Database_ptr_t &_database,
			const CommandLineOptions &_opts,
			const bool _checkTrueSolution
			);

	virtual ~InverseProblemSolver() {}

	/** Operator to solve this specific feasibility problem.
	 *
	 * @param _startingvalue start value or zero
	 * @return result container with success state and solutions
	 */
	GeneralMinimizer::ReturnValues operator()(
			const SpaceElement_ptr_t &_startingvalue)
	{ return operator()(_startingvalue, SpaceElement_ptr_t()); }

	/** Operator to solve this specific feasibility problem while
	 * comparing to a true solution.
	 *
	 * @param _startingvalue start value or zero
	 * @param _truesolution true solution for true error calculation
	 * @return result container with success state and solutions
	 */
	GeneralMinimizer::ReturnValues operator()(
			const SpaceElement_ptr_t &_startingvalue,
			const SpaceElement_ptr_t _truesolution);

	/** Returns specific name of the problem for output.
	 *
	 * @return name of the problem
	 */
	const std::string& getName() const
	{ return name; }

	/** Resets the feasibility problem.
	 *
	 */
	void clear();

	/** Finalizes the feasibility problem.
	 *
	 */
	void finish();

private:
	//!> ref to inverse problem
	InverseProblem_ptr_t inverseproblem;
	//!> database reference for storing iteration information
	Database_ptr_t database;
	//!> storing options on the manner of solving
	const CommandLineOptions &opts;
	//!> whether to calculate true solution in l2 case from SVD
	const bool checkTrueSolution;

private:
	//!> pre-constructed minimizer
	MinimizerFactory::instance_ptr_t minimizer;

	const std::string name;
};



#endif /* SOLVERS_INVERSEPROBLEMSOLVER_HPP_ */
