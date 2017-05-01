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

#include "Minimizations/Mappings/Specifics/DualRegularizedL1DualityMapping.hpp"
#include "Minimizations/types.hpp"

class SoftThresholdingMappingUnitTest;

/** This class implements the soft thresholding regularized l1 operator
 * where the regularization parameter \f$ lambda \f$ is kept fixed.
 *
 */
class SoftThresholdingMapping :
	public DualRegularizedL1DualityMapping
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

	/** Gets the one dual element of possibly many that minimizes the
	 * scalar product to the given other element.
	 *
	 * \param _x point for which to get dual element(s)
	 * \param _y other element to minimize scalar product against
	 * \return _Jx duality mapped \a _x for which \f$ < Jx, y > \f$ is minimal
	 */
	void getMinimumInfimum(
			const SpaceElement_ptr_t &_x,
			const SpaceElement_ptr_t &_y,
			SpaceElement_ptr_t &_Jx) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	const Mapping_ptr_t getAdjointMapping() const;
};


#endif /* SOFTTHRESHOLDINGMAPPING_HPP_ */
