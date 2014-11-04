/*
 * VectorProjection_BregmanDistanceToLine.cpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#include "VectorProjection_BregmanDistanceToLine.hpp"

#include <boost/log/trivial.hpp>
#include <cassert>

typedef typename MinimizationFunctional<Eigen::VectorXd>::array_type array_type;

BregmanDistanceToLine::BregmanDistanceToLine(
		const BregmanDistance &_distance,
		const Norm &_dualnorm,
		const PowerTypeDualityMapping &_J_q,
		const Eigen::VectorXd &_linevector,
		const Eigen::VectorXd &_tobeprojected,
		const double _powertype
		) :
		distance(_distance),
		dualnorm(_dualnorm),
		J_q(_J_q),
		linevector(_linevector),
		argvector(_tobeprojected),
		powertype(_powertype)
{
	// linevector normalized in respective norm
	assert( !linevector.isZero() );
	const_cast<Eigen::VectorXd &>(linevector) *= 1./dualnorm(_linevector);
}

double
BregmanDistanceToLine::function(
		const double &_value) const
{
	assert( fabs(dualnorm(linevector) - 1.) < std::numeric_limits<double>::epsilon()*1e2 );
	const Eigen::VectorXd Jy = _value*linevector;
	const double result = distance(argvector, Jy);
	return result;
}

const double
BregmanDistanceToLine::gradient(
		const double &_value) const
{
	assert( fabs(dualnorm(linevector) - 1.) < std::numeric_limits<double>::epsilon()*1e2 );
	const Eigen::VectorXd Jy = (_value)*linevector;
	const double grad = (J_q(Jy) - J_q(argvector)).transpose()
			* (linevector);
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
