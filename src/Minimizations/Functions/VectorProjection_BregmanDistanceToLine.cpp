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

typedef typename MinimizationFunctional< std::vector<double> >::array_type array_type;

BregmanDistanceToLine::BregmanDistanceToLine(
		const BregmanDistance &_distance,
		const Norm &_dualnorm,
		const PowerTypeDualityMapping &_J_q,
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
