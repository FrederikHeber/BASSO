/*
 * VectorProjection_BregmanDistanceToLine.cpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#include "VectorProjection_BregmanDistanceToLine.hpp"

#include <boost/log/trivial.hpp>


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
BregmanDistanceToLine::operator()(
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
BregmanDistanceToLine::convertToInternalType(
		double &_t,
		const gsl_vector * const x) const
{
	_t = gsl_vector_get(x, 0);
}

void
BregmanDistanceToLine::convertFromInternalType(
		const double &_t,
		gsl_vector * x) const
{
	gsl_vector_set(x, 0, _t);
}
