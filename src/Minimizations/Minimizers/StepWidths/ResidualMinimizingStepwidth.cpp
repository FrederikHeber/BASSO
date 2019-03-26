/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
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





