/*
 * InverseProblem.hpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */

#ifndef INVERSEPROBLEM_HPP_
#define INVERSEPROBLEM_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

/** This class represents all the required members of a given inverse
 * problem in an abstract form to the subsequent ..Minimizers. In effect,
 * it is a translator from the given parameters describing fully the inverse
 * problem to
 */
class InverseProblem
{
public:
	/** Constructor for class InverseProblem
	 *
	 * @param _A (Linear)Mapping to invert
	 * @param _y right-hand side to match
	 */
	InverseProblem(
			const Mapping_ptr_t &_A,
			const SpaceElement_ptr_t &_y
			);

	//!> mapping, defining problem from some space X to another space Y
	const Mapping_ptr_t A;
	//!> right-hand side of the problem in space Y
	const SpaceElement_ptr_t y;
	//!> sought-after solution in source space X
	SpaceElement_ptr_t x;
};



#endif /* INVERSEPROBLEM_HPP_ */
