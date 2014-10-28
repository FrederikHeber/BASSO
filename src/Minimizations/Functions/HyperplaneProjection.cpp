/*
 * HyperplaneProjection.cpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/log/trivial.hpp>

#include "HyperplaneProjection.hpp"

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
HyperplaneProjection::convertToInternalType(
		Eigen::VectorXd &_t,
		const gsl_vector * const x) const
{
	for (unsigned int i=0;i<_t.innerSize();++i)
		_t[i] = gsl_vector_get(x, i);
}

void
HyperplaneProjection::convertFromInternalType(
		const Eigen::VectorXd &_t,
		gsl_vector * x) const
{
	for (unsigned int i=0; i<_t.innerSize();++i)
		gsl_vector_set(x, i, _t[i]);
}

