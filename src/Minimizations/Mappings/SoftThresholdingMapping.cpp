/*
 * SoftThresholdingMapping.cpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#include "SoftThresholdingMapping.hpp"

#include "Math/Helpers.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

SoftThresholdingMapping::SoftThresholdingMapping(
		const NormedSpace_ptr_t &_NormedSpaceRef) :
	L1DualityMapping(_NormedSpaceRef, 2.),
	lambda(0.1)
{}

SoftThresholdingMapping::SoftThresholdingMapping(
		const NormedSpace_ptr_t &_NormedSpaceRef,
		const double _lambda) :
	L1DualityMapping(_NormedSpaceRef, 2.),
	lambda(_lambda)
{}

SpaceElement_ptr_t
SoftThresholdingMapping::operator()(
		const SpaceElement_ptr_t &_x) const
{
	SpaceElement_ptr_t result = _x->getSpace()->createElement();
	*result = operator()(_x->getVectorRepresentation());
	return result;
}

const Eigen::VectorXd
SoftThresholdingMapping::operator()(
		const Eigen::VectorXd &_x) const
{
	Eigen::VectorXd result(_x);
	for (unsigned int i=0;i<result.innerSize();++i)
		result[i] = fabs(result[i]) < lambda ?
				0. :
				(fabs(result[i])-lambda)*Helpers::sign(result[i]);
	return result;
}


