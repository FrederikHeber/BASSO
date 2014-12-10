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
#include "Minimizations/Mappings/SoftThresholdingMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"

/** This class implements a regularized l1 norm of the form
 * \f$ \lambda ||.||_1 + \tfrac 1 2 ||.||^2_2 \f$.
 *
 */
class DualRegularizedL1Norm : public LpNorm
{
public:
	/** Constructor for class Norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _ref reference to the space this norm is associated with
	 */
	DualRegularizedL1Norm(
			const NormedSpace_ptr_t& _ref,
			const double _lambda = 0.1) :
		LpNorm(_ref, 2.),
		softthresholder(_ref, _lambda)
	{}

	/** Evaluates the norm for a given \a _element.
	 *
	 * @param _element element of the space, whose norm to evaluated
	 * @return norm of \a element
	 */
	const double operator()(const SpaceElement_ptr_t &_x) const
	{
		assert( NormedSpaceRef == _x->getSpace() );
		const double value =
				LpNorm::operator()(softthresholder(_x));
		return value;
	}

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

private:
	//!> internal soft thresholding operator
	mutable SoftThresholdingMapping softthresholder;
};



#endif /* DUALREGULARIZEDL1NORM_HPP_ */
