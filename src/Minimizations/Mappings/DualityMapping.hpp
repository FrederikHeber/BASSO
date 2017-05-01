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

	/** Gets the one dual element of possibly many that minimizes the
	 * scalar product to the given other element.
	 *
	 * \param _x point for which to get dual element(s)
	 * \param _y other element to minimize scalar product against
	 * \return _Jx duality mapped \a _x for which \f$ < Jx, y > \f$ is minimal
	 */
	virtual void getMinimumInfimum(
			const SpaceElement_ptr_t &_x,
			const SpaceElement_ptr_t &_y,
			SpaceElement_ptr_t &_Jx) const = 0;
};


#endif /* DUALITYMAPPING_HPP_ */
