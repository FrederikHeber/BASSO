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
#include "Minimizations/DualityMapping.hpp"
#include "Minimizations/LpNorm.hpp"


double BregmanDistance::operator()(
		const Eigen::VectorXd &_x,
		const Eigen::VectorXd &_y,
		const double _power
		) const
{
	BOOST_LOG_TRIVIAL(trace)
			<< "Calculating Bregman distance between "
			<< _x.transpose() << " and " << _y.transpose();
	double result = 0.;
	if (p == LpNorm::Infinity) {
		result += _x.array().abs().maxCoeff();
		result += _y.array().abs().maxCoeff();
	} else {
		result += (1./q) * ::pow(lpnorm(_x), p);
		result += (1./p) * ::pow(lpnorm(_y), p);
	}
	result -= J_p(_x,_power).transpose() * _y;
	return result;
}
