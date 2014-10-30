/*
 * LInfinityDualityMapping.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef LINFINITYDUALITYMAPPING_HPP_
#define LINFINITYDUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"

class LInfinityDualityMappingUnitTest;

/** This class contains a duality mapping instance from a specific
 * l_infty space to its dual.
 *
 */
class LInfinityDualityMapping : public PowerTypeDualityMapping
{
	//!> grant unit test LInfinityDualityMappingUnitTest access to private parts.
	friend class LInfinityDualityMappingUnitTest;
public:
	/** Constructor for class LInfinityDualityMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _power power type of this duality mapping
	 */
	LInfinityDualityMapping(
			const NormedSpace_ptr_t &_NormedSpaceRef,
			const double _power) :
		PowerTypeDualityMapping(_NormedSpaceRef, _power)
	{}

	/** Evaluates duality mapping at \a _x.
	 *
	 * \param _x point where to evaluate
	 */
	virtual const Eigen::VectorXd operator()(
			const Eigen::VectorXd &_x) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	Mapping_ptr_t getAdjointMapping() const;
};

#endif /* LINFINITYDUALITYMAPPING_HPP_ */
