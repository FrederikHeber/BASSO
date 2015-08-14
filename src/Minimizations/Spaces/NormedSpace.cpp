/*
 * NormedSpace.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormedSpace.hpp"

#include <cassert>
#include <boost/bind.hpp>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/VectorSpaceOperationCounts.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

NormedSpace::NormedSpace(
		const unsigned int _dimension) :
	dimension(_dimension)
{}

NormedSpace::NormedSpace(
		const unsigned int _dimension,
		const Norm_ptr_t &_norm
		) :
	norm(_norm),
	dimension(_dimension)
{}

SpaceElement_ptr_t NormedSpace::createElement() const
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::ElementCreation>(
					opcounts.instance));
	SpaceElement_ptr_t newelement(new SpaceElement(Space));
	newelement->setSelfRef(newelement);
	return newelement;
}

const Mapping_ptr_t& NormedSpace::getDualityMapping() const
{
	assert (dualitymapping != NULL);
	return dualitymapping;
}
