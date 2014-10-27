/*
 * PowerTypeDualityMapping.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "PowerTypeDualityMapping.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

SpaceElement_ptr_t PowerTypeDualityMapping::operator()(
		const SpaceElement_ptr_t &_sourceelement,
		const double _power
		) const
{
	SpaceElement_ptr_t targetelement =
			_sourceelement->getSpace()->getDualSpace()->createElement();
	*targetelement = operator()(_sourceelement->getVectorRepresentation(), _power);
	return targetelement;
}
