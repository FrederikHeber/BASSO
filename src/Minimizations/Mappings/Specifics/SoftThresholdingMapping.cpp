/*
 * SoftThresholdingMapping.cpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#include "SoftThresholdingMapping.hpp"

#include "Math/Helpers.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/IllegalDualityMapping.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

SoftThresholdingMapping::SoftThresholdingMapping(
		const NormedSpace_weakptr_t &_NormedSpaceRef,
		const double _lambda) :
	DualRegularizedL1DualityMapping(_NormedSpaceRef, _lambda)
{}

void SoftThresholdingMapping::operator()(
		const SpaceElement_ptr_t &_x,
		SpaceElement_ptr_t &_Jx) const
{
	// start timing
	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	for (unsigned int i=0;i<_Jx->getSpace()->getDimension();++i)
		(*_Jx)[i] = fabs((*_x)[i]) < lambda ?
				0. :
				(fabs((*_x)[i])-lambda)*Helpers::sign((*_x)[i]);

	// finish timing
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	timing += timing_end - timing_start;
	++count;
}

const Mapping_ptr_t SoftThresholdingMapping::getAdjointMapping() const
{
	Mapping_ptr_t instance(
			new IllegalDualityMapping
	);
	return instance;
}

