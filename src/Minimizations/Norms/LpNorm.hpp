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
	LpNorm(const NormedSpace_weakptr_t& _ref,
			const double _p) :
		Norm(_ref),
		p(_p),
		p_as_int(p)
	{
		if ((p <= 1.) || (p == std::numeric_limits<double>::infinity()))
			throw NormIllegalValue_exception()
				<< NormIllegalValue_name("p");

		if (fabs((double)p_as_int - p) > BASSOTOLERANCE)
			const_cast<int &>(p_as_int) = 0;
	}
	~LpNorm() {}

	/** Getter for the p value of a possible lp norm.
	 *
	 * @return p value: 0 - not an lp norm, else - p of lp norm
	 */
	virtual const double getPvalue() const
	{ return p; }

	bool isSmooth() const
	{ return true; }

protected:

	/** Evaluates the norm for a given \a _element.
	 *
	 * @param _element element of the space, whose norm to evaluated
	 * @return norm of \a element
	 */
	const double internal_operator(const SpaceElement_ptr_t &_x) const
	{
		assert( getSpace() == _x->getSpace() );
		double value = 0.;
		if (p_as_int == 0)
			for (unsigned int i=0;i<_x->getSpace()->getDimension();++i)
				value += ::pow(fabs((*_x)[i]), p);
		else
			for (unsigned int i=0;i<_x->getSpace()->getDimension();++i)
				value += ::pow(fabs((*_x)[i]), p_as_int);
		return ::pow(value, 1./p);
	}

private:
	//!> p value for norm
	const double p;

	//!> states whether p is actually integer value (or is 0 otherwise)
	const int p_as_int;
};


#endif /* LPNORM_HPP_ */
