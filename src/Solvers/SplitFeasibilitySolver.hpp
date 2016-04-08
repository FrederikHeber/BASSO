/*
 * SplitFeasibilitySolver.hpp
 *
 *  Created on: Nov 19, 2015
 *      Author: heber
 */

#ifndef SOLVERS_SPLITFEASIBILITYSOLVER_HPP_
#define SOLVERS_SPLITFEASIBILITYSOLVER_HPP_

#include "BassoConfig.h"

#include <deque>

#include "Minimizations/types.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints.hpp"
#include "Solvers/FeasibilityProblem.hpp"

class CommandLineOptions;

/** The SplitFeasibilitySolver is akin to the InverseProblemSolver but taking
 * auxiliary constraints into account.
 *
 * It only differs by taking additional auxiliary constraints into account
 * and solving the inverse problem and each constraint in turn until
 * convergence.
 */
class SplitFeasibilitySolver : public FeasibilityProblem
{
public:
	/** Cstor for class SplitFeasibilitySolver.
	 *
	 * @param _opts options defining manner of inverse problem solving
	 */
	SplitFeasibilitySolver(
			const CommandLineOptions &_opts);

	virtual ~SplitFeasibilitySolver() {};

	/** Register a FeasibilityProblem into the qeue of this
	 * solver.
	 *
	 * @param _fp Feasibility problem
	 */
	void registerFeasibilityProblem(
			FeasibilityProblem::ptr_t &_fp);

	/** Register auxiliary constraints to fulfill.
	 *
	 * @param _auxiliary_constraints auxiliary constraints
	 */
	void registerAuxiliaryConstraints(
			const AuxiliaryConstraints::ptr_t &_auxiliary_constraints
			);

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
	 * Here, we don't do anything. clear() and finish() are used
	 * inside operator() already.
	 *
	 */
	void clear()
	{}

	/** Finalizes the feasibility problem.
	 *
	 * Here, we don't do anything. clear() and finish() are used
	 * inside operator() already.
	 *
	 */
	void finish()
	{}

private:
	//!> storing options on the manner of solving
	const CommandLineOptions &opts;

	//!> typedef for deque with problems
	typedef std::deque<FeasibilityProblem::ptr_t> problems_t;
	//!> problems to be project onto solution manifold cyclically
	problems_t problems;

	const std::string name;
};



#endif /* SOLVERS_SPLITFEASIBILITYSOLVER_HPP_ */
