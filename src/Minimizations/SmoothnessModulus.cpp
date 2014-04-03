/*
 * SmoothnessModulus.cpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#include "SmoothnessModulus.hpp"

//#include <boost/log/trivial.hpp>
#include <cmath>

#include "MinimizationExceptions.hpp"

SmoothnessModulus::SmoothnessModulus(const double _p) :
	p(_p)
{
	if (p <= 1.) throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("p");
}

double SmoothnessModulus::operator()(const double _value)
{
//	BOOST_LOG_TRIVIAL(trace) <<
//			"Calculating modulus of smoothness at " << _value;

	if ((_value <= 0.) && (_value > 1.)) throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("_value");

	double result = 0.;
	if ((p > 1.) && (p <= 2.)) {
		const double value_p = ::pow(_value, p);
		result = ::pow( 1. + value_p, 1./p) - 1.;
	} else if (p > 2.) {
		result = .5*(p-1.)*_value*_value;
	}

	return result;
}
