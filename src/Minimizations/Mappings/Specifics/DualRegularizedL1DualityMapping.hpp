/*
 * DualRegularizedL1DualityMapping.hpp
 *
 *  Created on: May 02, 2017
 *      Author: heber
 */

#ifndef DUALREGULARIZEDL1DUALITYMAPPING_HPP_
#define DUALREGULARIZEDL1DUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>

#include "Minimizations/Mappings/Specifics/L1DualityMapping.hpp"
#include "Minimizations/types.hpp"

/** This class defines the interface for all regularized l1 duality
 * mappings that require a lambda regularization parameter, according to
 * [Schoepfer 2012] which is required for the dual norm and duality
 * mapping to the (squared) regularized l1 norm.
 */
class DualRegularizedL1DualityMapping :
	public L1DualityMapping
{
public:
	/** Default constructor of class DualRegularizedL1DualityMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _lambda soft thresholding parameter
	 */
	DualRegularizedL1DualityMapping(
			const NormedSpace_weakptr_t &_NormedSpaceRef,
			const double _lambda = 0.1) :
		L1DualityMapping(_NormedSpaceRef, 2.),
		lambda(_lambda),
		count(0),
		timing(boost::chrono::nanoseconds(0))
	{}

	//!> expose overloaded operator method from base class
	using Mapping::operator();

	/** Evaluates for the given \a _x the soft thresholding result with respect
	 * to \a _lambda.
	 * @param _x vector to relative shrink
	 * @param _Jx componentwise relative shrunk of \a _x by \a lambda
	 */
	virtual void operator()(
			const SpaceElement_ptr_t &_x,
			SpaceElement_ptr_t &_Jx) const = 0;

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

protected:
	//!> soft thresholding parameter
	const double lambda;

	//!> number of times the operator was called
	mutable unsigned int count;

	//!> total runtime spent on executing this operator
	mutable boost::chrono::nanoseconds timing;
};


#endif /* DUALREGULARIZEDL1DUALITYMAPPING_HPP_ */
