/*
 * IllegalDualityMapping.cpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#include "IllegalDualityMapping.hpp"

#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

IllegalDualityMapping::IllegalDualityMapping() :
	LpDualityMapping(NormedSpaceFactory::getDummySpace(), 2.)
{}

const NormedSpace_ptr_t IllegalDualityMapping::getSourceSpace() const
{
	// we just throw as this function must now be called
	throw MinimizationIllegalValue_exception()
		<< MinimizationIllegalValue_name("illegally called");
}

const NormedSpace_ptr_t IllegalDualityMapping::getTargetSpace() const
{
	// we just throw as this function must now be called
	throw MinimizationIllegalValue_exception()
		<< MinimizationIllegalValue_name("illegally called");
}

void IllegalDualityMapping::operator()(
		const SpaceElement_ptr_t &_sourceelement,
		SpaceElement_ptr_t &_destelement
		) const
{
	// we just throw as this function must now be called
	throw MinimizationIllegalValue_exception()
		<< MinimizationIllegalValue_name("illegally called");
}

const Mapping_ptr_t IllegalDualityMapping::getAdjointMapping() const
{
	// we just throw as this function must now be called
	throw MinimizationIllegalValue_exception()
		<< MinimizationIllegalValue_name("illegally called");
}


