/*
 * Specifics/SoftThresholdingMapping.hpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#ifndef SOFTTHRESHOLDINGMAPPING_HPP_
#define SOFTTHRESHOLDINGMAPPING_HPP_

#include "BassoConfig.h"

#include "Minimizations/Mappings/Specifics/L1DualityMapping.hpp"
#include "Minimizations/types.hpp"

class SoftThresholdingMappingUnitTest;

class SoftThresholdingMapping :
	public L1DualityMapping
{
	//!> grant unit test access to private parts
	friend class SoftThresholdingMappingUnitTest;
public:
	/** Default constructor of class SoftThresholdingMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 */
	SoftThresholdingMapping(const NormedSpace_ptr_t &_NormedSpaceRef);

	/** Default constructor of class SoftThresholdingMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _lambda soft thresholding parameter
	 */
	SoftThresholdingMapping(
			const NormedSpace_ptr_t &_NormedSpaceRef,
			const double _lambda);

	/** Evaluates for the given \a _x the soft thresholding result with respect
	 * to \a _lambda.
	 * @param _x vector to soft-threshold
	 * @return componentwise soft threshold of \a _x by \a lambda
	 */
	const SpaceElement_ptr_t operator()(
			const SpaceElement_ptr_t &_x) const;

	/** Setter for soft thresholding parameter \a lambda.
	 *
	 * @param _lambda new value for parameter
	 */
	void setLambda(const double _lambda) const
	{
		const_cast<double &>(lambda) = _lambda;
	}

	/** Getter for the regularization parameter.
	 *
	 * @return regularization parameter
	 */
	const double getLambda() const
	{ return lambda; }

private:
	//!> soft thresholding parameter
	const double lambda;
};


#endif /* SOFTTHRESHOLDINGMAPPING_HPP_ */
