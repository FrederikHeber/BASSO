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
#include "Minimizations/Mappings/DualityMapping.hpp"
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
	double resultSKP = 0.;
	// check whether space is smooth
	if (_x->getSpace()->getNorm()->isSmooth()) {
		resultSKP =_xdual * _y;
	} else {
		// if not smooth, duality mapping is not single-value and hence we have
		// to evaluate the infimum over all possible dual elements
		const DualityMapping &mapping =
				static_cast<DualityMapping &>(*_x->getSpace()->getDualityMapping());
		SpaceElement_ptr_t xdual_mininf = _xdual->getSpace()->createElement();
		mapping.getMinimumInfimum(_x, _y, xdual_mininf);
		resultSKP = xdual_mininf * _y;
	}
	const double resultX = (1./Helpers::ConjugateValue(power)) * ::pow(norm(_x), power);
	const double resultY = (1./power) * ::pow(norm(_y), power);
	double result = resultX+resultY-resultSKP;
	LOG(debug, "Calculating Bregman distance between " << _x << " and " << _y
			<< " as " << resultX << "+" << resultY << "-" << resultSKP << "="
			<< result);
	if (_x->getSpace()->getNorm()->isSmooth()) {
		resultSKP =_xdual * _y;
		double old_result = resultX+resultY-resultSKP;
		LOG(debug, "Calculating former Bregman distance between " << _x << " and " << _y
				<< " as " << resultX << "+" << resultY << "-" << resultSKP << "="
				<< old_result);
	}
	return result;
}
