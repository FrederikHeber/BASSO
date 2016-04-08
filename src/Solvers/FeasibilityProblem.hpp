/*
 * FeasibilityProblem.hpp
 *
 *  Created on: Nov 20, 2015
 *      Author: heber
 */

#ifndef SOLVERS_FEASIBILITYPROBLEM_HPP_
#define SOLVERS_FEASIBILITYPROBLEM_HPP_

#include "BassoConfig.h"

#include <string>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"

/** This describe the interface for a single part of a split
 * feasibility problem.
 */
struct FeasibilityProblem
{
	//!> shared ptr typedef to safely wrap a problem instance
	typedef boost::shared_ptr<FeasibilityProblem> ptr_t;

	/** Operator to solve this specific feasibility problem.
	 *
	 * @param _startingvalue start value or zero
	 * @return result container with success state and solutions
	 */
	virtual GeneralMinimizer::ReturnValues operator()(
			const SpaceElement_ptr_t &_startingvalue) = 0;

	/** Operator to solve this specific feasibility problem while
	 * comparing to a true solution.
	 *
	 * @param _startingvalue start value or zero
	 * @param _truesolution true solution for true error calculation
	 * @return result container with success state and solutions
	 */
	virtual GeneralMinimizer::ReturnValues operator()(
			const SpaceElement_ptr_t &_startingvalue,
			const SpaceElement_ptr_t _truesolution) = 0;

	/** Returns specific name of the problem for output.
	 *
	 * @return name of the problem
	 */
	virtual const std::string& getName() const = 0;

	/** Resets the feasibility problem.
	 *
	 * This is used for SplitFeasibilityProblems where the
	 * problem remains the same and only the starting values
	 * changes upon the multiple runs
	 *
	 */
	virtual void clear() = 0;

	/** Finalizes the feasibility problem.
	 *
	 * This is used for SplitFeasibilityProblems where the
	 * problem remains the same and only the starting values
	 * changes upon the multiple runs. In the end we need to
	 * accumulate some values from the database.
	 *
	 */
	virtual void finish() = 0;
};

#endif /* SOLVERS_FEASIBILITYPROBLEM_HPP_ */
