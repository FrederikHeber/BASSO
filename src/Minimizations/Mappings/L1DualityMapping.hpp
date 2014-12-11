/*
 * L1DualityMapping.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef L1DUALITYMAPPING_HPP_
#define L1DUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"

class L1DualityMappingUnitTest;

/** This class contains a duality mapping instance from a specific
 * l_1 space to its dual.
 *
 */
class L1DualityMapping : public PowerTypeDualityMapping
{
	//!> grant unit test L1DualityMappingUnitTest access to private parts.
	friend class L1DualityMappingUnitTest;
public:
	/** Constructor for class L1DualityMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _power power type of this duality mapping
	 */
	L1DualityMapping(
			const NormedSpace_ptr_t &_NormedSpaceRef,
			const double _power) :
		PowerTypeDualityMapping(_NormedSpaceRef, _power)
	{}

	/** Evaluates duality mapping at \a _x.
	 *
	 * \param _x point where to evaluate
	 */
	virtual const SpaceElement_ptr_t operator()(
			const SpaceElement_ptr_t &_x) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	const Mapping_ptr_t getAdjointMapping() const;
};

#endif /* L1DUALITYMAPPING_HPP_ */
