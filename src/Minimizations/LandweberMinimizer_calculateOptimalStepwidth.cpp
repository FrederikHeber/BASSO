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

#include "MinimizationExceptions.hpp"
#include "Functions/ResidualFunctional.hpp"

double LandweberMinimizer::calculateOptimalStepwidth(
		const Eigen::VectorXd &_x,
		const Eigen::VectorXd &_dualx,
		const Eigen::VectorXd &_u,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y,
		const double _alpha) const
{
	double alpha = _alpha;
	ResidualFunctional res(
			_x,	// x_n
			_dualx,	// x^*_n
			_u, // u_n
			_A, // A
			_y, // y
			*this // landweber
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



