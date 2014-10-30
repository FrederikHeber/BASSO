/*
 * LandweberMinimizer_calculateOptimalStepwidth.cpp
 *
 *  Created on: Apr 14, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LandweberMinimizer.hpp"

#include <Eigen/Dense>
#include <boost/bind.hpp>
#include <boost/math/tools/minima.hpp>
#include <limits>

#include "Minimizations/Functions/ResidualFunctional.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"

double LandweberMinimizer::calculateOptimalStepwidth(
		const InverseProblem_ptr_t &_problem,
		const SpaceElement_ptr_t &_dualx,
		const SpaceElement_ptr_t &_u,
		const double _alpha) const
{
	double alpha = _alpha;
	ResidualFunctional res(
			_problem,
			_dualx,	// x^*_n
			_u, // u_n
			*this // minimizer
			);
	double minval = -1000.;
	double maxval = 1000.;
	int bits = 32;
	boost::uintmax_t maxiter = 100;
	std::pair<double, double> minpair =
			boost::math::tools::brent_find_minima(
					boost::bind(
							&ResidualFunctional::operator(),
							boost::cref(res),
							_1),
					minval,
					maxval,
					bits,
					maxiter);
	alpha = minpair.first;

	return alpha;
}



