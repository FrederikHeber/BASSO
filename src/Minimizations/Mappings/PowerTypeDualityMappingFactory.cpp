/*
 * PowerTypeDualityMappingFactory.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "PowerTypeDualityMappingFactory.hpp"

#include <limits>

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/L1DualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LpDualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LInfinityDualityMapping.hpp"

const Mapping_ptr_t
PowerTypeDualityMappingFactory::createInstance(
		const NormedSpace_weakptr_t &_NormedSpaceRef,
		const double _power)
{
	PowerTypeDualityMapping *mapping = NULL;
	const NormedSpace_ptr_t Space(_NormedSpaceRef);
	const double p =
			Space->getNorm()->getPvalue();
	if (p == std::numeric_limits<double>::infinity()) {
		mapping = new LInfinityDualityMapping(Space,_power);
	} else if (p > 1.) {
		mapping = new LpDualityMapping(Space,_power);
	} else if (p == 1.) {
		mapping = new L1DualityMapping(Space,_power);
	}
	return Mapping_ptr_t(mapping);
}
