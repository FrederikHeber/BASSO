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
 * VectorProjection_BregmanDistanceToLine.cpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#include "VectorProjection_BregmanDistanceToLine.hpp"

#include <boost/log/trivial.hpp>
#include <cassert>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Norms/Norm.hpp"

typedef typename MinimizationFunctional< std::vector<double> >::array_type array_type;

BregmanDistanceToLine::BregmanDistanceToLine(
		const BregmanDistance &_distance,
		const Norm &_dualnorm,
		const Mapping &_J_q,
		const SpaceElement_ptr_t &_linevector,
		const SpaceElement_ptr_t &_tobeprojected,
		const double _powertype
		) :
		distance(_distance),
		dualnorm(_dualnorm),
		J_q(_J_q),
		linevector(_linevector),
		argvector(_tobeprojected),
		powertype(_powertype),
		linevector_norm(dualnorm(linevector))
{}

double
BregmanDistanceToLine::function(
		const double &_value) const
{
	const SpaceElement_ptr_t Jy = (1./linevector_norm)*_value*linevector;
	const double result = distance(argvector, Jy);
	return result;
}

const double
BregmanDistanceToLine::gradient(
		const double &_value) const
{
	const SpaceElement_ptr_t Jy = (1./linevector_norm)*(_value)*linevector;
	const double grad = (J_q(Jy) - J_q(argvector))
			* (linevector)/linevector_norm;
	return grad;
}

void
BregmanDistanceToLine::convertInternalTypeToArrayType(
		const double &_t,
		array_type & _x
		) const
{
	assert( _x.size() == 1);
	_x[0] = _t;
}

void
BregmanDistanceToLine::convertArrayTypeToInternalType(
		const array_type & _x,
		double &_t
		) const
{
	assert( _x.size() == 1);
	_t = _x[0];
}
