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
	L1Norm(
			const NormedSpace_ptr_t& _ref) :
		Norm(_ref)
	{}

	const double operator()(const SpaceElement_ptr_t &_x) const
	{
		assert( NormedSpaceRef == _x->getSpace() );
		double value = 0.;
		for (unsigned int i=0;i<_x->getSpace()->getDimension();++i)
			value += fabs((*_x)[i]);
		return value;
	}

	/** Getter for the p value of a possible lp norm.
	 *
	 * @return p value: 0 - not an lp norm, else - p of lp norm
	 */
	virtual const double getPvalue() const
	{ return 1.; }
};


#endif /* L1NORM_HPP_ */
