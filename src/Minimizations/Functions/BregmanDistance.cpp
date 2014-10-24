/*
 * BregmanDistance.cpp
 *
 *  Created on: Apr 25, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "BregmanDistance.hpp"

#include <cmath>
#include <Eigen/Dense>

#include "Log/Logging.hpp"
#include "Minimizations/DualityMappings/DualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"

double BregmanDistance::operator()(
		const Eigen::VectorXd &_x,
		const Eigen::VectorXd &_y,
		const double _power
		) const
{
	const Eigen::VectorXd dual_x = J_p(_x,_power).transpose();
	return operator()(_x,_y, dual_x, _power);
}

double BregmanDistance::operator()(
		const Eigen::VectorXd &_x,
		const Eigen::VectorXd &_y,
		const Eigen::VectorXd &_xdual,
		const double _power
		) const
{

	if ((_power != 0.) && (_power <= 1.))
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("_power");
	BOOST_LOG_TRIVIAL(trace)
			<< "Calculating Bregman distance between "
			<< _x.transpose() << " and " << _y.transpose();
	double result = 0.;
	if (p == LpNorm::Infinity) {
		result += _x.array().abs().maxCoeff();
		result += _y.array().abs().maxCoeff();
	} else {
		result += ((_power-1.)/_power) * ::pow(lpnorm(_x), _power);
		result += (1./_power) * ::pow(lpnorm(_y), _power);
	}
	result -= ScalarVectorProduct(_xdual, _y);
	return result;
}
