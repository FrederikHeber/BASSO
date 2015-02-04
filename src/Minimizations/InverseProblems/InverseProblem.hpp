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
 * problem to the elements required for a solving algorithm.
 */
class InverseProblem
{
public:
	/** Constructor for class InverseProblem
	 *
	 * \note We store source and target space for both convenience and
	 * because we need to keep the strong shared_ptr's somewhere as
	 * SpaceElement and LinearMapping have only weak_ptr's to their
	 * associated spaces.
	 *
	 * @param _A (Linear)Mapping to invert
	 * @param _SourceSpace ptr to source space
	 * @param _TargetSpace ptr to targetspace
	 * @param _y right-hand side to match
	 */
	InverseProblem(
			const Mapping_ptr_t &_A,
			const NormedSpace_ptr_t &_SourceSpace,
			const NormedSpace_ptr_t &_TargetSpace,
			const SpaceElement_ptr_t &_y
			);

	//!> mapping, defining problem from some space X to another space Y
	const Mapping_ptr_t A;
	//!> mapping, defining problem from some space X to another space Y
	const Mapping_ptr_t A_t;
	//!> sought-for solution in space X
	SpaceElement_ptr_t x;
	//!> right-hand side of the problem in space Y
	const SpaceElement_ptr_t y;
	//!> source space
	const NormedSpace_ptr_t SourceSpace;
	//!> dual source space
	const NormedSpace_ptr_t DualSourceSpace;
	//!> target space
	const NormedSpace_ptr_t TargetSpace;
	//!> dual target space
	const NormedSpace_ptr_t DualTargetSpace;
};



#endif /* INVERSEPROBLEM_HPP_ */
