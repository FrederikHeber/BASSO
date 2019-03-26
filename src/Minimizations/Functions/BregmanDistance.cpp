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
 * BregmanDistance.cpp
 *
 *  Created on: Apr 25, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "BregmanDistance.hpp"

#include <cmath>
#include <limits>

#include "Log/Logging.hpp"
#include "Math/Helpers.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Norms/Norm.hpp"

BregmanDistance::BregmanDistance(
		const Norm &_norm,
		const Mapping &_J_p,
		const double _power) :
			power(_power),
			norm(_norm),
			J_p(_J_p)
{
	if ((power != 0.) && (power <= 1.))
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("power");
}

BregmanDistance::BregmanDistance(
		const InverseProblem_ptr_t &_problem) :
			power(_problem->SourceSpace->getDualityMapping()->getPower()),
			norm(*_problem->SourceSpace->getNorm()),
			J_p(*_problem->SourceSpace->getDualityMapping())
{
	if ((power != 0.) && (power <= 1.))
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("power");
}

double BregmanDistance::operator()(
		const SpaceElement_ptr_t &_x,
		const SpaceElement_ptr_t &_y
		) const
{
	const SpaceElement_ptr_t dual_x = J_p(_x);
	return operator()(_x,_y, dual_x);
}

double BregmanDistance::operator()(
		const SpaceElement_ptr_t &_x,
		const SpaceElement_ptr_t &_y,
		const SpaceElement_ptr_t &_xdual
		) const
{

//	LOG(trace, "Calculating Bregman distance between " << _x << " and " << _y);
	double result = 0.;
	result += (1./Helpers::ConjugateValue(power)) * ::pow(norm(_x), power);
	result += (1./power) * ::pow(norm(_y), power);
	result -= _xdual * _y;
	return result;
}
