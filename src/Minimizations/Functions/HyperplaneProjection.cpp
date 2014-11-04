/*
 * HyperplaneProjection.cpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/log/trivial.hpp>

#include "HyperplaneProjection.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"

typedef typename MinimizationFunctional<Eigen::VectorXd>::array_type array_type;

HyperplaneProjection::HyperplaneProjection(
	BregmanProjectionFunctional &_bregman,
	const Eigen::VectorXd &_x,
	const Eigen::MatrixXd &_U,
	const Eigen::VectorXd &_alpha
	) :
		bregman(_bregman),
		x(_x),
		U(_U),
		alpha(_alpha)
{}

HyperplaneProjection::HyperplaneProjection(
	BregmanProjectionFunctional &_bregman,
	const SpaceElement_ptr_t &_x,
	const Eigen::MatrixXd &_U,
	const Eigen::VectorXd &_alpha
	) :
		bregman(_bregman),
		x(_x->getVectorRepresentation()),
		U(_U),
		alpha(_alpha)
{}

double
HyperplaneProjection::operator()(
		const Eigen::VectorXd &_value) const
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

const Eigen::VectorXd
HyperplaneProjection::gradient(
		const Eigen::VectorXd &_value) const
{
	const Eigen::VectorXd grad =
			bregman.gradient(
					_value,
					x,
					U,
					alpha);
	BOOST_LOG_TRIVIAL(trace)
		<< "grad() evaluates to " << grad.transpose();
	return grad;
}

void
HyperplaneProjection::convertToArrayType(
		const Eigen::VectorXd &_t,
		array_type & _x
		) const
{
	for (unsigned int i=0; i<_t.innerSize();++i)
		_x[i] = _t[i];
}

void
HyperplaneProjection::convertFromArrayType(
		const array_type & _x,
		Eigen::VectorXd &_t
		) const
{
	for (unsigned int i=0;i<_t.innerSize();++i)
		_t[i] = _x[i];
}

