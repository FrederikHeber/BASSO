/*
 * L1Norm.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef L1NORM_HPP_
#define L1NORM_HPP_

#include "BassoConfig.h"

#include <cassert>
#include <cmath>

#include "Norm.hpp"

#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

/** This class implements the l_1 norm.
 *
 */
class L1Norm : public Norm
{
public:
	/** Constructor of class L1Norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 */
	L1Norm(const NormedSpace_weakptr_t& _ref) :
		Norm(_ref)
	{}

	/** Getter for the p value of a possible lp norm.
	 *
	 * @return p value: 0 - not an lp norm, else - p of lp norm
	 */
	virtual const double getPvalue() const
	{ return 1.; }

	/** We return smooth although the l_1 norm is not smooth but a
	 * suitable selection exists.
	 *
	 * @return true
	 */
	bool isSmooth() const
	{ return true; }

protected:

	const double internal_operator(const SpaceElement_ptr_t &_x) const
	{
		assert( getSpace() == _x->getSpace() );
		const Eigen::VectorXd &vector = RepresentationAdvocate::get(_x);
		return vector.array().abs().sum();
	}
};


#endif /* L1NORM_HPP_ */
