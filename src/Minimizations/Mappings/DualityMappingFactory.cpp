/*
 * DualityMappingFactory.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */


#include <Minimizations/Mappings/DualityMappingFactory.hpp>
#include "BassoConfig.h"

#include <limits>

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/L1DualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LpDualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LInfinityDualityMapping.hpp"

const Mapping_ptr_t
DualityMappingFactory::createInstance(
		const NormedSpace_weakptr_t _NormedSpaceRef,
		const double _power)
{
	PowerTypeDualityMapping *mapping = NULL;
	const double p =
			NormedSpace_ptr_t(_NormedSpaceRef)->getNorm()->getPvalue();
	if (p == std::numeric_limits<double>::infinity()) {
		mapping = new LInfinityDualityMapping(_NormedSpaceRef,_power);
	} else if (p > 1.) {
		mapping = new LpDualityMapping(_NormedSpaceRef,_power);
	} else if (p == 1.) {
		mapping = new L1DualityMapping(_NormedSpaceRef,_power);
	}
	return Mapping_ptr_t(mapping);
}
