/*
 * DualRegularizedL1Norm.hpp
 *
 *  Created on: Oct 31, 2014
 *      Author: heber
 */

#ifndef DUALREGULARIZEDL1NORM_HPP_
#define DUALREGULARIZEDL1NORM_HPP_

#include "BassoConfig.h"

#include <cassert>
#include "Minimizations/Mappings/Specifics/RelativeShrinkageMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

/** This class implements the dual of the regularized l1 norm of the form
 * \f$ \lambda ||.||_1 + \tfrac 1 2 ||.||^2_2 \f$ which is
 * \f$ || S_{\lambda} (.) ||_2 \f$ .
 *
 */
class DualRegularizedL1Norm : public Norm
{
public:
	/** Constructor for class Norm.
	 *
	 * \note The internal RelativeShrinkageMapping brings the given \a _x
	 * into its dual space. Hence, we set the internal L2-norm
	 * onto this space as well.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _ref reference to the space this norm is associated with
	 */
	DualRegularizedL1Norm(
			const NormedSpace_weakptr_t& _ref,
			const double _lambda = 0.1) :
		Norm(_ref),
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
		const double value = l2norm(softthresholder(_x));
		return value;
	}

private:
	//!> internal soft thresholding operator
	mutable RelativeShrinkageMapping softthresholder;

	//!> internal l2 norm
	const LpNorm l2norm;
};



#endif /* DUALREGULARIZEDL1NORM_HPP_ */
