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

//	BOOST_LOG_TRIVIAL(trace)
//			<< "Calculating Bregman distance between "
//			<< _x << " and " << _y;
	double result = 0.;
	result += (1./Helpers::ConjugateValue(power)) * ::pow(norm(_x), power);
	result += (1./power) * ::pow(norm(_y), power);
	result -= _xdual * _y;
	return result;
}
