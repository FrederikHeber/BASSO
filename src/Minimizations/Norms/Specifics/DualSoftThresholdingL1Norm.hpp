/*
 * DualSoftThresholdingL1Norm.hpp
 *
 *  Created on: Oct 31, 2014
 *      Author: heber
 */

#ifndef DUALSOFTTHRESHOLDINGL1NORM_HPP_
#define DUALSOFTTHRESHOLDINGL1NORM_HPP_

#include "BassoConfig.h"

#include <cassert>
#include "Minimizations/Mappings/Specifics/SoftThresholdingMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/Specifics/DualRegularizedL1Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

/** This class implements the dual of the regularized l1 norm of the form
 * \f$ \sqrt{ \frac 1 2 ||.||^2_1 + \frac \lambda 2 ||.||^2_2 } \f$ which is
 * \f$ || \sqrt{ c_{\lambda}^2(.) + \frac 1 {\lambda} || S_{\lambda} (.) ||^2_2 } \f$ .
 *
 * see [Schoepfer, 2012].
 */
class DualSoftThresholdingL1Norm : public DualRegularizedL1Norm
{
public:
	/** Constructor for class Norm.
	 *
	 * \note The internal RelativeShrinkageMapping brings the given \a _x
	 * into its dual space. Hence, we set the internal L2-norm
	 * onto this space as well.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _lambda regularization parameter
	 */
	DualSoftThresholdingL1Norm(
			const NormedSpace_weakptr_t& _ref,
			const double _lambda = 0.1) :
		DualRegularizedL1Norm(_ref),
		softthresholder(_ref, _lambda),
		l2norm(softthresholder.getTargetSpace(), 2.)
	{}

	/** Setter for soft thresholding parameter \a lambda.
	 *
	 * @param _lambda new value for parameter
	 */
	void setLambda(const double _lambda) const
	{
		softthresholder.setLambda(_lambda);
	}

	/** Getter for the regularization parameter.
	 *
	 * @return regularization parameter
	 */
	const double getLambda() const
	{ return softthresholder.getLambda(); }

protected:

	/** Evaluates the norm for a given \a _x.
	 *
	 *  @param _x element of the space, whose norm to evaluated
	 * @return norm of \a _x
	 */
	const double internal_operator(const SpaceElement_ptr_t &_x) const
	{
		assert( getSpace() == _x->getSpace() );
		double value = 0.;
		value += ::pow(getLambda(), 2);
		value += ::pow(l2norm(softthresholder(_x)), 2)/getLambda();
		return sqrt(value);
	}

private:
	//!> internal soft thresholding operator
	mutable SoftThresholdingMapping softthresholder;

	//!> internal l2 norm
	const LpNorm l2norm;
};



#endif /* DUALSOFTTHRESHOLDINGL1NORM_HPP_ */
