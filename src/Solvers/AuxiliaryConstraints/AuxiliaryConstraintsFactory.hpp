/*
 * AuxiliaryConstraintsFactory.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef SOLVERS_AUXILIARYCONSTRAINTS_AUXILIARYCONSTRAINTSFACTORY_HPP_
#define SOLVERS_AUXILIARYCONSTRAINTS_AUXILIARYCONSTRAINTSFACTORY_HPP_

#include "BassoConfig.h"

#include <map>

#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints.hpp"

/** Constructs stopping criteria for given string to be parsed by internal
 * lexer or by elemental type.
 *
 */
struct AuxiliaryConstraintsFactory
{
	/** Constructor for AuxiliaryConstraintsFactory.
	 *
	 */
	AuxiliaryConstraintsFactory();

	//!> enumeration of all stopping criteria
	enum AuxiliaryConstraints_types
	{
		Nonnegative,	//!< components must be non-negative
		Nonpositive,	//!< components must be non-positive
		Unity,			//!< components must have norm of 1
		MAX_AuxiliaryConstraints_types
	};

	//!> enumeration of all combination types
	enum Combination_types
	{
		NoCombination=0,
		Combination_AND,
		MAX_Combination_types
	};

	/** Creates a combined stopping criterion out of a given string in
	 * \a _criteria_line.
	 *
	 * @param _criteria_line string to convert to stopping criteria
	 * @return combined stopping criterion
	 */
	AuxiliaryConstraints::ptr_t create(
			const std::string &_criteria_line);

	/** Creates a single basic stopping criterion.
	 *
	 * @param _type desired type for the stopping criterion
	 * @return type
	 */
	AuxiliaryConstraints::ptr_t createCriterion(
			const enum AuxiliaryConstraints_types _type);

	/** Creates a single basic combination of stopping criteria.
	 *
	 * @param _type desired type for the stopping criterion
	 * @param _left first stopping criterion
	 * @param _right second stopping criterion (set NULL for NOT)
	 * @return combined criterion
	 */
	AuxiliaryConstraints::ptr_t
	createCombination(
			const enum Combination_types _type,
			const AuxiliaryConstraints::ptr_t &_left,
			const AuxiliaryConstraints::ptr_t &_right);

	/** Converts stopping criterion name to type.
	 *
	 * Returns MAX_AuxiliaryConstraints_types when name not found in map
	 *
	 * @param _name name of stopping criterion
	 * @return type of stopping criterion
	 */
	enum AuxiliaryConstraints_types getCriterionTypeForName(const std::string &_name);

	/** Converts combination name to type.
	 *
	 * Returns MAX_Combination_types when name not found in map
	 *
	 * @param _name name of combination
	 * @return type of combination
	 */
	enum Combination_types getCombinationTypeForName(const std::string &_name);

private:
	//!> typedef for stopping criterion name to type map
	typedef std::map<std::string, AuxiliaryConstraints_types> StringToCriterionTypeMap_t;
	//!> stopping criterion name to type map
	StringToCriterionTypeMap_t StringToCriterionTypeMap;
	//!> typedef for combiner name to type map
	typedef std::map<std::string, Combination_types> StringToCombinationTypeMap_t;
	//!> combiner name to type map
	StringToCombinationTypeMap_t StringToCombinationTypeMap;
};


#endif /* SOLVERS_AUXILIARYCONSTRAINTS_AUXILIARYCONSTRAINTSFACTORY_HPP_ */
