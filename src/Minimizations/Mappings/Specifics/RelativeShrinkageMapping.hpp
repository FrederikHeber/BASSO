/*
 * Specifics/RelativeShrinkageMapping.hpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#ifndef SOFTTHRESHOLDINGMAPPING_HPP_
#define SOFTTHRESHOLDINGMAPPING_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>

#include "Minimizations/Mappings/Specifics/L1DualityMapping.hpp"
#include "Minimizations/types.hpp"

class RelativeShrinkageMappingUnitTest;

/** This class implements a Relative Shrinkage operator according to
 * [Schoepfer 2012] which is required for the dual norm and duality
 * mapping to the (squared) regularized l1 norm.
 */
class RelativeShrinkageMapping :
	public L1DualityMapping
{
	//!> grant unit test access to private parts
	friend class RelativeShrinkageMappingUnitTest;
public:
	/** Default constructor of class RelativeShrinkageMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 */
	RelativeShrinkageMapping(const NormedSpace_weakptr_t &_NormedSpaceRef);

	virtual ~RelativeShrinkageMapping() {}

	/** Default constructor of class RelativeShrinkageMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _lambda soft thresholding parameter
	 */
	RelativeShrinkageMapping(
			const NormedSpace_weakptr_t &_NormedSpaceRef,
			const double _lambda);

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

	/** Returns the number of times the operator() was called.
	 *
	 * @return number of calls
	 */
	const unsigned int getCount() const
	{ return count; }

	/** Returns the total runtime the program spent so far on
	 * its operator().
	 *
	 * @return runtime summed over all calls
	 */
	const boost::chrono::nanoseconds getTiming() const
	{ return timing; }

	/** Calculates the relative shrinkage coefficient \f$ c_{\lambda} (x) \f$
	 * depending on the given \a x.
	 *
	 * @param _x element to calculate shrinkage for
	 * @return relative shrinkage coefficient
	 */
	const double getRelativeShrinkage(const SpaceElement_ptr_t &_x) const;

protected:
	//!> soft thresholding parameter
	const double lambda;

	//!> number of times the operator was called
	mutable unsigned int count;

	//!> total runtime spent on executing this operator
	mutable boost::chrono::nanoseconds timing;
};


#endif /* SOFTTHRESHOLDINGMAPPING_HPP_ */
