/*
 * Specifics/LpDualityMapping.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#ifndef LPDUALITYMAPPING_HPP_
#define LPDUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"

class LpDualityMappingUnitTest;

/** This class contains a duality mapping instance from a specific
 * l_p space with 1 < p < infty to its dual.
 *
 */
class LpDualityMapping : public PowerTypeDualityMapping
{
	//!> grant unit test LpDualityMappingUnitTest access to private parts.
	friend class LpDualityMappingUnitTest;
public:
	/** Constructor for class LpDualityMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _power power type of this duality mapping
	 */
	LpDualityMapping(
			const NormedSpace_ptr_t &_NormedSpaceRef,
			const double _power);

	/** Mapping function.
	 *
	 * @param _sourceelement element to map/transform
	 * @return new transformed/mapped element
	 */
	const SpaceElement_ptr_t operator()(
			const SpaceElement_ptr_t &_sourceelement
			) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	const Mapping_ptr_t getAdjointMapping() const;
};

#endif /* LPDUALITYMAPPING_HPP_ */
