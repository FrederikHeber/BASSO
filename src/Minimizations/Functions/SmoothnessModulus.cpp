/*
 * SmoothnessModulus.cpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "SmoothnessModulus.hpp"

//#include <boost/log/trivial.hpp>
#include <cmath>

#include "Minimizations/Minimizers/MinimizationExceptions.hpp"

SmoothnessModulus::SmoothnessModulus(const double _p) :
	p(_p)
{
	if (p <= 1.) throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("p");
}

double SmoothnessModulus::operator()(const double _value) const
{
//	LOG(trace,
//			"Calculating modulus of smoothness at " << _value);

	if ((_value <= 0.) && (_value > 1.)) throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("_value");

	double result = 0.;
	if ((p > 1.) && (p <= 2.)) {
		const double value_p = ::pow(_value, p);
		result = value_p/p;
				// ::pow( 1. + value_p, 1./p) - 1.;
	} else if (p > 2.) {
		result = .5*(p-1.)*_value*_value;
//				::pow(
//				.5*(::pow(1. + _value, p) + ::pow(fabs(1. - _value), p)),
//				1./p) - 1.; //
	}

	return result;
}
