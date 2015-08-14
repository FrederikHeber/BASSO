/*
 * ResidualMinimizingStepwidth.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "ResidualMinimizingStepwidth.hpp"

#include <boost/bind.hpp>
#include <boost/math/tools/minima.hpp>
#include <limits>

#include "Minimizations/Functions/ResidualFunctional.hpp"

const double ResidualMinimizingStepwidth::operator()(
		const SpaceElement_ptr_t &_dualx,
		const SpaceElement_ptr_t &_u,
		const SpaceElement_ptr_t &_solution,
		const SpaceElement_ptr_t &_residual,
		const double _residuum,
		const double _TolX,
		const double _alpha
		) const
{
	double alpha = _alpha;
	ResidualFunctional res(
			problem,
			_dualx,	// x^*_n
			_u, // u_n
			residualizer // residual calculator
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





