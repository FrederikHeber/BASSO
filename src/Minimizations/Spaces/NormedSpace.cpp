/*
 * NormedSpace.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormedSpace.hpp"

#include <boost/bind.hpp>
#include <Eigen/Dense>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/VectorSpaceOperationCounts.hpp"

NormedSpace::NormedSpace(
		const unsigned int _dimension) :
	dimension(_dimension)
{}

NormedSpace::NormedSpace(
		const unsigned int _dimension,
		const Norm_ptr_t &_norm,
		const constructDualityMapping_t &_constructDualityMapping
		) :
	norm(_norm),
	constructDualityMapping(_constructDualityMapping),
	dimension(_dimension)
{}

SpaceElement_ptr_t NormedSpace::createElement() const
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::ElementCreation>(
					opcounts.instance));
	SpaceElement_ptr_t newelement(new SpaceElement(getSpace()));
	newelement->setSelfRef(newelement);
	return newelement;
}

const Mapping_ptr_t& NormedSpace::getDualityMapping() const
{
	if (dualitymapping == NULL) {
		const_cast<Mapping_ptr_t &>(dualitymapping) =
				constructDualityMapping();
	}
	return dualitymapping;
}

