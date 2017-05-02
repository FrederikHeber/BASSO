/*
 * RelativeShrinkageMapping.hpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#ifndef RELATIVESHRINKAGEMAPPING_HPP_
#define RELATIVESHRINKAGEMAPPING_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>

#include "Minimizations/Mappings/Specifics/DualRegularizedL1DualityMapping.hpp"
#include "Minimizations/types.hpp"

class RelativeShrinkageMappingUnitTest;

/** This class implements a Relative Shrinkage operator according to
 * [Schoepfer 2012] which is required for the dual norm and duality
 * mapping to the (squared) regularized l1 norm.
 */
class RelativeShrinkageMapping :
	public DualRegularizedL1DualityMapping
{
	//!> grant unit test access to private parts
	friend class RelativeShrinkageMappingUnitTest;
public:
	/** Default constructor of class RelativeShrinkageMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _lambda soft thresholding parameter
	 */
	RelativeShrinkageMapping(
			const NormedSpace_weakptr_t &_NormedSpaceRef,
			const double _lambda = 0.1);

	//!> expose overloaded operator method from base class
	using Mapping::operator();

	/** Evaluates for the given \a _x the soft thresholding result with respect
	 * to \a _lambda.
	 * @param _x vector to relative shrink
	 * @param _Jx componentwise relative shrunk of \a _x by \a lambda
	 */
	void operator()(
			const SpaceElement_ptr_t &_x,
			SpaceElement_ptr_t &_Jx) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	const Mapping_ptr_t getAdjointMapping() const;

	/** Calculates the relative shrinkage coefficient \f$ c_{\lambda} (x) \f$
	 * depending on the given \a x.
	 *
	 * @param _x element to calculate shrinkage for
	 * @return relative shrinkage coefficient
	 */
	const double getRelativeShrinkage(const SpaceElement_ptr_t &_x) const;
};


#endif /* RELATIVESHRINKAGEMAPPING_HPP_ */
