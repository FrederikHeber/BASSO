/*
 * LInfinityNorm.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef LINFINITYNORM_HPP_
#define LINFINITYNORM_HPP_

#include "BassoConfig.h"

#include <limits>

#include "Norm.hpp"

#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

/** This class implements the l_infinity norm.
 *
 */
class LInfinityNorm : public Norm
{
public:
	/** Constructor of class LInfinityNorm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 */
	LInfinityNorm(const NormedSpace_weakptr_t& _ref) :
		Norm(_ref)
	{}

	/** Getter for the p value of a possible lp norm.
	 *
	 * @return p value: 0 - not an lp norm, else - p of lp norm
	 */
	virtual const double getPvalue() const
	{ return std::numeric_limits<double>::infinity(); }

protected:

	/** Evaluates the norm for a given \a _element.
	 *
	 * @param _element element of the space, whose norm to evaluated
	 * @return norm of \a element
	 */
	const double internal_operator(const SpaceElement_ptr_t &_x) const
	{
		// infinity norm
		assert( getSpace() == _x->getSpace() );
		return _x->getMaxCoefficientAndIndex().first;
	}
};


#endif /* LINFINITYNORM_HPP_ */
