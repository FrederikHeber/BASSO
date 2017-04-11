/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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
