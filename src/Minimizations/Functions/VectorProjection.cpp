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
 * VectorProjection.cpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "VectorProjection.hpp"

#include <utility>

#include "Log/Logging.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizerFactory.hpp"
#include "Minimizations/Functions/VectorProjection_BregmanDistanceToLine.hpp"

VectorProjection::VectorProjection(
		const Norm &_lpnorm,
		const Mapping &_J_p,
		const double _p
		) :
		lpnorm(_lpnorm),
		J_p(_J_p),
		p(_p)
{}

const std::pair<double, double>
VectorProjection::operator()(
		const SpaceElement_ptr_t &_projectedonto,
		const SpaceElement_ptr_t &_tobeprojected,
		const double _Tol) const
{
	BregmanDistance distance(
			lpnorm,
			J_p,
			p
			);


	const unsigned int dim = 1;
	double tmin = 0.;
	double result = 0.;
	{

		// tmin=fminunc(@(t) BregmanFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
		const BregmanDistanceToLine distancefunctional(
				distance,
				lpnorm,
				J_p,
				_projectedonto,
				_tobeprojected,
				p);
		FunctionalMinimizer<double>::ptr_t functionminimizer =
				FunctionalMinimizerFactory::create<double>(
						dim,
						distancefunctional,
						tmin);

//		const unsigned int inner_iterations =
				(*functionminimizer)(dim, _Tol, tmin);
		result = functionminimizer->getCurrentOptimumValue();

		LOG(trace, "tmin is " << tmin);
	}

	return std::make_pair( result, tmin);
}
