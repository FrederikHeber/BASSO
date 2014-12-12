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

typedef VectorSpaceOperationCounts::TimeKeeper TimeKeeper;

NormedSpace::NormedSpace(
		const unsigned int _dimension) :
	dimension(_dimension)
{}

NormedSpace::NormedSpace(
		const unsigned int _dimension,
		const Norm_ptr_t &_norm,
		const Mapping_ptr_t &_dualitymapping
		) :
	norm(_norm),
	dualitymapping(_dualitymapping),
	dimension(_dimension)
{}

SpaceElement_ptr_t NormedSpace::createElement() const
{
	TimeKeeper(opcounts.ElementCreation);
	SpaceElement_ptr_t newelement(new SpaceElement(getSpace()));
	newelement->setSelfRef(newelement);
	return newelement;
}
