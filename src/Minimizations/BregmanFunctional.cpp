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
		const LpNorm &_lpdualnorm,
		const DualityMapping &_J_q
		) :
	lpdualnorm(_lpdualnorm),
	J_q(_J_q)
{}

double BregmanFunctional::operator()(
		const Eigen::VectorXd &_t,
		const Eigen::VectorXd &_dualx,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha,
		const double _q
		)
{
	// x=x-U*t;
	const Eigen::VectorXd resx = _dualx - _U * _t;
	// fval=1/q*norm(x,p)^q+alpha'*t;
	const double fval =
			1./(double)_q * ::pow(lpdualnorm(resx), _q)
			+ _alpha.transpose() * _t;
	return fval;
}

Eigen::VectorXd BregmanFunctional::gradient(
		const Eigen::VectorXd &_t,
		const Eigen::VectorXd &_dualx,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha,
		const double _q
		)
{
	const Eigen::VectorXd resx = _dualx - _U * _t;
	const Eigen::VectorXd gval =
			_alpha -
			_U.transpose() * J_q(resx, _q);

	return gval;
}
