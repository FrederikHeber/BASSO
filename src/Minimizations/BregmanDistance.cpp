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

#include "Minimizations/DualityMapping.hpp"
#include "Minimizations/LpNorm.hpp"


double BregmanDistance::operator()(
		const Eigen::VectorXd &_x,
		const Eigen::VectorXd &_y
		) const
{
	double result = 0.;
	if (p == LpNorm::Infinity) {
		result += _x.array().abs().maxCoeff();
		result += _y.array().abs().maxCoeff();
	} else {
		result += (1./q) * ::pow(lpnorm(_x), p);
		result += (1./p) * ::pow(lpnorm(_y), p);
	}
	result -= J_p(_x,p).transpose() * _y;
	return result;
}
