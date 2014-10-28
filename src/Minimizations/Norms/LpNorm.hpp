/*
 * LpNorm.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef LPNORM_HPP_
#define LPNORM_HPP_

#include "BassoConfig.h"

#include <cmath>
#include <Eigen/Dense>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/NormExceptions.hpp"

class LpNorm : public Norm
{
public:
	LpNorm(const double _p) :
		p(_p)
	{
		if (p < 0.)
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
//		assert( NormedSpaceRef == x->getSpace() );
		return operator()(_x->getVectorRepresentation());
	}

	const double operator()(const Eigen::VectorXd &_x) const
	{
		if (p > 0. ) {
			double value = 0.;
			for (unsigned int i=0;i<_x.innerSize();++i)
				value += ::pow(fabs(_x[i]), p);
			return ::pow(value, 1./p);
		} else {
			// infinity norm
			return _x.lpNorm<Eigen::Infinity>();
		}
	}

	enum { Infinity = 0};

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
