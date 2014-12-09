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

const SpaceElement_ptr_t
SoftThresholdingMapping::operator()(
		const SpaceElement_ptr_t &_x) const
{
	SpaceElement_ptr_t result = getTargetSpace()->createElement();
	for (unsigned int i=0;i<result->getSpace()->getDimension();++i)
		(*result)[i] = fabs((*_x)[i]) < lambda ?
				0. :
				(fabs((*_x)[i])-lambda)*Helpers::sign((*_x)[i]);
	return result;
}


