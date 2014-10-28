/*
 * SoftThresholdingOperator.cpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#include "SoftThresholdingOperator.hpp"

#include "Math/Helpers.hpp"

SoftThresholdingOperator::SoftThresholdingOperator() :
	L1DualityMapping(2.),
	lambda(0.1)
{}

SoftThresholdingOperator::SoftThresholdingOperator(
		const double _lambda) :
	L1DualityMapping(2.),
	lambda(_lambda)
{}

const Eigen::VectorXd
SoftThresholdingOperator::operator()(
		const Eigen::VectorXd &_x) const
{
	Eigen::VectorXd result(_x);
	for (unsigned int i=0;i<result.innerSize();++i)
		result[i] = fabs(result[i]) < lambda ?
				0. :
				(fabs(result[i])-lambda)*Helpers::sign(result[i]);
	return result;
}


