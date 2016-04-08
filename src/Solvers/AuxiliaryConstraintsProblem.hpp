/*
 * AuxiliaryConstraintsProblem.hpp
 *
 *  Created on: Nov 20, 2015
 *      Author: heber
 */

#ifndef SOLVERS_AUXILIARYCONSTRAINTSPROBLEM_HPP_
#define SOLVERS_AUXILIARYCONSTRAINTSPROBLEM_HPP_

#include "BassoConfig.h"

#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints.hpp"
#include "Solvers/FeasibilityProblem.hpp"

/** This class is an adaptor for a given auxiliary constraints
 * to be called as a FeasibilityProblem.
 */
class AuxiliaryConstraintsProblem : public FeasibilityProblem
{
public:
	/** Cstor of AuxiliaryConstraintsProblem.
	 *
	 */
	AuxiliaryConstraintsProblem(
			const AuxiliaryConstraints::ptr_t &_ac) :
		ac(_ac),
		name("AuxiliaryConstraints")
	{}

	virtual ~AuxiliaryConstraintsProblem() {};

	/** Operator to solve this specific feasibility problem.
	 *
	 * @param _startingvalue start value or zero
	 * @return result container with success state and solutions
	 */
	GeneralMinimizer::ReturnValues operator()(
			const SpaceElement_ptr_t &_startingvalue)
	{ return work(_startingvalue); }

	/** Operator to solve this specific feasibility problem while
	 * comparing to a true solution.
	 *
	 * @param _startingvalue start value or zero
	 * @param _truesolution true solution for true error calculation
	 * @return result container with success state and solutions
	 */
	GeneralMinimizer::ReturnValues operator()(
			const SpaceElement_ptr_t &_startingvalue,
			const SpaceElement_ptr_t _truesolution)
	{ return work(_startingvalue); }

	/** Returns specific name of the problem for output.
	 *
	 * @return name of the problem
	 */
	const std::string& getName() const
	{ return name; }

	/** Resets the feasibility problem.
	 *
	 */
	void clear()
	{}

	/** Finalizes the feasibility problem.
	 *
	 */
	void finish()
	{}

private:
	/** Helper function cause I don't know how to properly handle overloaded
	 * virtual functions with default arguments.
	 *
	 * @param _startingvalue start value or zero
	 * @return result container with success state and solutions
	 */
	GeneralMinimizer::ReturnValues work(
			const SpaceElement_ptr_t &_startingvalue)
	{
		GeneralMinimizer::ReturnValues result;

		result.m_solution = _startingvalue->getSpace()->createElement();
		*result.m_solution = _startingvalue;
		result.m_dual_solution =
				_startingvalue->getSpace()->getDualSpace()->createElement();
		result.m_dual_solution->setZero();
		(*ac)(result.m_solution);
		result.residuum = (_startingvalue - result.m_solution)->Norm();

		return result;
	}

private:
	//!> auxiliary constraints which this problem fulfills.
	AuxiliaryConstraints::ptr_t ac;

	const std::string name;
};



#endif /* SOLVERS_AUXILIARYCONSTRAINTSPROBLEM_HPP_ */
