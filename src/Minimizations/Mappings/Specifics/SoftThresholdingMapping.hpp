/*
 * Specifics/SoftThresholdingMapping.hpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#ifndef SOFTTHRESHOLDINGMAPPING_HPP_
#define SOFTTHRESHOLDINGMAPPING_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>

#include "Minimizations/Mappings/Specifics/RegularizedL1DualityMapping.hpp"
#include "Minimizations/types.hpp"

class SoftThresholdingMappingUnitTest;

/** This class implements the soft thresholding regularized l1 operator
 * where the regularization parameter \f$ lambda \f$ is kept fixed.
 *
 */
class SoftThresholdingMapping :
	public RegularizedL1DualityMapping
{
	//!> grant unit test access to private parts
	friend class SoftThresholdingMappingUnitTest;
public:
	/** Default constructor of class SoftThresholdingMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _lambda soft thresholding parameter
	 */
	SoftThresholdingMapping(
			const NormedSpace_weakptr_t &_NormedSpaceRef,
			const double _lambda = 0.1);

	//!> expose overloaded operator method from base class
	using Mapping::operator();

	/** Evaluates for the given \a _x the soft thresholding result with respect
	 * to \a _lambda.
	 * @param _x vector to soft-threshold
	 * @param _Jx componentwise soft threshold of \a _x by \a lambda
	 */
	void operator()(
			const SpaceElement_ptr_t &_x,
			SpaceElement_ptr_t &_Jx) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	const Mapping_ptr_t getAdjointMapping() const;
};


#endif /* SOFTTHRESHOLDINGMAPPING_HPP_ */
