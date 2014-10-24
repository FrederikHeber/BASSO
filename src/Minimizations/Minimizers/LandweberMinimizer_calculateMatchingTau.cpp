/*
 * LandweberMinimizer_calculateMatchingTau.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LandweberMinimizer.hpp"

#include <boost/bind.hpp>
#include <boost/log/trivial.hpp>
#include <boost/math/tools/minima.hpp>
#include <limits>

#include "Minimizations/Functions/SmoothnessFunctional.hpp"
#include "Minimizations/MinimizationExceptions.hpp"


double LandweberMinimizer::calculateMatchingTau(
		const SmoothnessModulus &_modul,
		const double _lambda
		) const
{
	SmoothnessFunctional smoothness(_modul, _lambda);
	double tau = .5;
	double minval = -1000.;
	double maxval = 1000.;
	int bits = 32;
	boost::uintmax_t maxiter = 100;
	std::pair<double, double> minpair =
			boost::math::tools::brent_find_minima(
					boost::bind(
							&SmoothnessFunctional::operator(),
							boost::cref(smoothness),
							_1),
					minval,
					maxval,
					bits,
					maxiter);
	tau = minpair.first;

	BOOST_LOG_TRIVIAL(trace)
		<< "Matching tau from modulus of smoothness is " << tau;
	BOOST_LOG_TRIVIAL(trace)
		<< "Counter-check: rho(tau)/tau = "
		<< _modul(tau)/tau << ", lambda = " << _lambda;

	return tau;
}


