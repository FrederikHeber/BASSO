/*
 * StoppingCriteriaFactory.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERIAFACTORY_HPP_
#define MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERIAFACTORY_HPP_

#include "BassoConfig.h"

#include <map>

#include "Minimizations/Minimizers/StoppingCriteria/StoppingArguments.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion.hpp"

/** Constructs stopping criteria for given string to be parsed by internal
 * lexer or by elemental type.
 *
 */
struct StoppingCriteriaFactory
{
	/** Constructor for StoppingCriteriaFactory.
	 *
	 */
	StoppingCriteriaFactory();

	//!> enumeration of all stopping criteria
	enum StoppingCriteria_types
	{
		IterationCount,
		RelativeResiduum,
		Residuum,
		Walltime,
		MAX_StoppingCriteria_types
	};

	//!> enumeration of all combination types
	enum Combination_types
	{
		NoCombination=0,
		Combination_AND,
		Combination_OR,
		Combination_NOT,
		MAX_Combination_types
	};

	/** Creates a combined stopping criterion out of a given string in
	 * \a _criteria_line.
	 *
	 * @param _criteria_line string to convert to stopping criteria
	 * @param _args arguments for the stopping criterion
	 * @return combined stopping criterion
	 */
	StoppingCriterion::ptr_t create(
			const std::string &_criteria_line,
			const StoppingArguments &_args);

	/** Creates a single basic stopping criterion.
	 *
	 * @param _type desired type for the stopping criterion
	 * @param _args arguments for the stopping criterion
	 * @return type
	 */
	StoppingCriterion::ptr_t createCriterion(
			const enum StoppingCriteria_types _type,
			const StoppingArguments &_args);

	/** Creates a single basic combination of stopping criteria.
	 *
	 * @param _type desired type for the stopping criterion
	 * @param _left first stopping criterion
	 * @param _right second stopping criterion (set NULL for NOT)
	 * @return combined criterion
	 */
	StoppingCriterion::ptr_t
	createCombination(
			const enum Combination_types _type,
			const StoppingCriterion::ptr_t &_left,
			const StoppingCriterion::ptr_t &_right);

	/** Converts stopping criterion name to type.
	 *
	 * Returns MAX_StoppingCriteria_types when name not found in map
	 *
	 * @param _name name of stopping criterion
	 * @return type of stopping criterion
	 */
	enum StoppingCriteria_types getCriterionTypeForName(const std::string &_name);

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
	typedef std::map<std::string, StoppingCriteria_types> StringToCriterionTypeMap_t;
	//!> stopping criterion name to type map
	StringToCriterionTypeMap_t StringToCriterionTypeMap;
	//!> typedef for combiner name to type map
	typedef std::map<std::string, Combination_types> StringToCombinationTypeMap_t;
	//!> combiner name to type map
	StringToCombinationTypeMap_t StringToCombinationTypeMap;
};


#endif /* MINIMIZATIONS_MINIMIZERS_STOPPINGCRITERIA_STOPPINGCRITERIAFACTORY_HPP_ */
