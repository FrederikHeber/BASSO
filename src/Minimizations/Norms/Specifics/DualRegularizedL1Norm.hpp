/*
 * DualRegularizedL1Norm.hpp
 *
 *  Created on: Apr 25, 2017
 *      Author: heber
 */

#ifndef DUALREGULARIZEDL1NORM_HPP_
#define DUALREGULARIZEDL1NORM_HPP_

#include "BassoConfig.h"

#include <cassert>
#include "Minimizations/Mappings/Specifics/RelativeShrinkageMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/Specifics/RelativeShrinkageCoefficient.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

/** This class defines the interface for the dual of the regularized l1 norm
 * of the form \f$ \sqrt{ \frac 1 2 ||.||^2_1 + \frac \lambda 2 ||.||^2_2 } \f$
 * which is \f$ || \sqrt{ c_{\lambda}^2(.) + \frac 1 {\lambda} || S_{\lambda} (.) ||^2_2 } \f$ .
 *
 * see [Schoepfer, 2012].
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
	 */
	DualRegularizedL1Norm(
			const NormedSpace_weakptr_t& _ref) :
		Norm(_ref)
	{}

	/** Setter for soft thresholding parameter \a lambda.
	 *
	 * @param _lambda new value for parameter
	 */
	virtual void setLambda(const double _lambda) const = 0;

	/** Getter for the regularization parameter.
	 *
	 *
	 * @return regularization parameter
	 */
	virtual const double getLambda() const = 0;

	bool isSmooth() const
	{ return true; }

protected:

	/** Evaluates the norm for a given \a _x.
	 *
	 *  @param _x element of the space, whose norm to evaluated
	 * @return norm of \a _x
	 */
	virtual const double internal_operator(const SpaceElement_ptr_t &_x) const = 0;
};



#endif /* DUALREGULARIZEDL1NORM_HPP_ */
