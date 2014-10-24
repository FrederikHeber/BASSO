/*
 * SmoothnessModulus.hpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#ifndef SMOOTHNESSMODULUS_HPP_
#define SMOOTHNESSMODULUS_HPP_

#include "BassoConfig.h"

/** Functor for the modulus of smoothness in Lp spaces.
 *
 * \note For the moment we just implemented the asymptotic formulas given in
 * [Schöpfer et al., 2006]
 *
 */
class SmoothnessModulus
{
public:
	SmoothnessModulus(const double _p);
	~SmoothnessModulus() {}

	/** Evaluates modulus at \a _value.
	 *
	 * \param _value
	 * \return rho of \a _value
	 */
	double operator()(const double _value) const;

private:
	//!> p value for norm
	const double p;
};


#endif /* SMOOTHNESSMODULUS_HPP_ */
