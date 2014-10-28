/*
 * RegularizedL1Norm.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef REGULARIZEDL1NORM_HPP_
#define REGULARIZEDL1NORM_HPP_

#include "BassoConfig.h"

#include <cassert>
#include "Minimizations/Norms/L1Norm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"

/** This class implements a regularized l1 norm of the form
 * \f$ \lambda ||.||_1 + \tfrac 1 2 ||.||^2_2 \f$.
 *
 */
class RegularizedL1Norm : public L1Norm
{
public:
	/** Constructor for class Norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 */
	RegularizedL1Norm(
//			const boost::shared_ptr<NormedSpace> _ref,
			const double _lambda = 0.1) :
		l2norm(2.),
		lambda(_lambda)
	{}

	/** Evaluates the norm for a given \a _element.
	 *
	 * @param _element element of the space, whose norm to evaluated
	 * @return norm of \a element
	 */
	const double operator()(const SpaceElement_ptr_t &_x) const
	{
//		assert( NormedSpaceRef == x->getSpace() );
		return operator()(_x->getVectorRepresentation());
	}

	const double operator()(const Eigen::VectorXd &_x) const
	{
		double value = 0.;
		value += lambda * L1Norm::operator()(_x);
		value += .5 * ::pow(l2norm(_x),2);
		return value;
	}

private:
	//!> internal l2 norm
	LpNorm l2norm;
	//!> regularization parameter lambda
	const double lambda;
};



#endif /* REGULARIZEDL1NORM_HPP_ */