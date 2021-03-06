/*
 * DualityMapping.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef DUALITYMAPPING_HPP_
#define DUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

/** This class defines mappings from a space to its dual space with
 * elements being the gradient of the norm of the source space.
 *
 * Hence, only one space needs to be given as the other is automatically
 * the space's dual space.
 */
class DualityMapping : public Mapping
{
public:
	/** Constructor for class DualityMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 */
	DualityMapping(
			const NormedSpace_weakptr_t &_NormedSpaceRef
			) :
		Mapping(
				_NormedSpaceRef,
				NormedSpace_ptr_t(_NormedSpaceRef)->getDualSpace())
	{}

	virtual ~DualityMapping() {}

	//!> expose overloaded operator method from base class
	using Mapping::operator();

	/** Evaluates duality mapping at \a _x.
	 *
	 * \param _x point where to evaluate
	 * \param _Jx duality mapped \a _x
	 */
	virtual void operator()(
			const SpaceElement_ptr_t &_x,
			SpaceElement_ptr_t &_Jx) const = 0;
};


#endif /* DUALITYMAPPING_HPP_ */
