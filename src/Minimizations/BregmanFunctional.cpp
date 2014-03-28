/*
 * BregmanFunctional.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "BregmanFunctional.hpp"

#include <cmath>
#include <Eigen/Dense>

#include "Minimizations/DualityMapping.hpp"
#include "Minimizations/LpNorm.hpp"

BregmanFunctional::BregmanFunctional(
		const double _p,
		const double _tolerance
		) :
	p(_p),
	lpnorm(p),
	J_p(p)
{
	J_p.setTolerance(_tolerance);
}

double BregmanFunctional::operator()(
		const Eigen::VectorXd &_t,
		const Eigen::VectorXd &_x,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha,
		const unsigned int _q
		)
{
	// x=x-U*t;
	const Eigen::VectorXd resx = _x - _U * _t;
	// fval=1/q*norm(x,p)^q+alpha'*t;
	const double fval =
			1./(double)_q * ::pow(lpnorm(resx), _q)
			+ _alpha.transpose() * _t;
	return fval;
}

Eigen::VectorXd BregmanFunctional::gradient(
		const Eigen::VectorXd &_t,
		const Eigen::VectorXd &_x,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha,
		const unsigned int _q
		)
{
	const Eigen::VectorXd resx = _x - _U * _t;
	const Eigen::VectorXd gval =
			_alpha -
			_U.transpose() * J_p(resx, _q);

	return gval;
}
