/*
 * HyperplaneProjection.cpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "HyperplaneProjection.hpp"

#include "Log/Logging.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"

typedef typename MinimizationFunctional< std::vector<double> >::array_type array_type;

HyperplaneProjection::HyperplaneProjection(
	BregmanProjectionFunctional &_bregman,
	const SpaceElement_ptr_t &_x,
	const std::vector<SpaceElement_ptr_t> &_U,
	const std::vector<double> &_alpha
	) :
		bregman(_bregman),
		x(_x),
		U(_U),
		alpha(_alpha)
{}

double
HyperplaneProjection::function(
		const std::vector<double> &_value) const
{
	const double returnvalue =
			bregman(
					_value,
					x,
					U,
					alpha);
	BOOST_LOG_TRIVIAL(trace)
		<< "func() evaluates to " << returnvalue;
	return returnvalue;
}

const std::vector<double>
HyperplaneProjection::gradient(
		const std::vector<double> &_value) const
{
	const std::vector<double> grad =
			bregman.gradient(
					_value,
					x,
					U,
					alpha);
	std::stringstream outputstream;
	outputstream << grad;
	BOOST_LOG_TRIVIAL(trace)
		<< "grad() evaluates to " << outputstream.str();
	return grad;
}

void
HyperplaneProjection::convertInternalTypeToArrayType(
		const std::vector<double> &_t,
		array_type & _x
		) const
{
	_x = _t;
}

void
HyperplaneProjection::convertArrayTypeToInternalType(
		const array_type & _x,
		std::vector<double> &_t
		) const
{
	_t = _x;
}

