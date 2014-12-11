/*
 * LpNorm.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef LPNORM_HPP_
#define LPNORM_HPP_

#include "BassoConfig.h"

#include <cassert>
#include <cmath>

#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/NormExceptions.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

class LpNorm : public Norm
{
public:
	/** Constructor of class L1Norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _p p value of the norm
	 */
	LpNorm(
			const NormedSpace_ptr_t& _ref,
			const double _p) :
		Norm(_ref),
		p(_p)
	{
		if ((p <= 1.) || (p == std::numeric_limits<double>::infinity()))
			throw NormIllegalValue_exception()
				<< NormIllegalValue_name("p");

	}
	~LpNorm() {}

	/** Evaluates the norm for a given \a _element.
	 *
	 * @param _element element of the space, whose norm to evaluated
	 * @return norm of \a element
	 */
	const double operator()(const SpaceElement_ptr_t &_x) const
	{
		assert( NormedSpaceRef == _x->getSpace() );
		double value = 0.;
		for (unsigned int i=0;i<_x->getSpace()->getDimension();++i)
			value += ::pow(fabs((*_x)[i]), p);
		return ::pow(value, 1./p);
	}

	/** Getter for the p value of a possible lp norm.
	 *
	 * @return p value: 0 - not an lp norm, else - p of lp norm
	 */
	virtual const double getPvalue() const
	{ return p; }

private:
	//!> p value for norm
	const double p;
};


#endif /* LPNORM_HPP_ */
