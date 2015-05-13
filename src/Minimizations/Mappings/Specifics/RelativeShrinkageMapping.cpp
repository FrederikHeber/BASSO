/*
 * RelativeShrinkageMapping.cpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#include "RelativeShrinkageMapping.hpp"

#include "Math/Helpers.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/IllegalDualityMapping.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

RelativeShrinkageMapping::RelativeShrinkageMapping(
		const NormedSpace_weakptr_t &_NormedSpaceRef) :
	L1DualityMapping(_NormedSpaceRef, 2.),
	lambda(0.1),
	count(0),
	timing(boost::chrono::nanoseconds(0))
{}

RelativeShrinkageMapping::RelativeShrinkageMapping(
		const NormedSpace_weakptr_t &_NormedSpaceRef,
		const double _lambda) :
	L1DualityMapping(_NormedSpaceRef, 2.),
	lambda(_lambda),
	count(0),
	timing(boost::chrono::nanoseconds(0))
{}

const SpaceElement_ptr_t
RelativeShrinkageMapping::operator()(
		const SpaceElement_ptr_t &_x) const
{
	// start timing
	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	SpaceElement_ptr_t result = getTargetSpace()->createElement();
	for (unsigned int i=0;i<result->getSpace()->getDimension();++i)
		(*result)[i] = fabs((*_x)[i]) < lambda ?
				0. :
				(fabs((*_x)[i])-lambda)*Helpers::sign((*_x)[i]);

	// finish timing
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	timing += timing_end - timing_start;
	++count;

	return result;
}

const Mapping_ptr_t RelativeShrinkageMapping::getAdjointMapping() const
{
	Mapping_ptr_t instance(
			new IllegalDualityMapping
	);
	return instance;
}

