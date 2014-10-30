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
#include "Minimizations/Mappings/L1DualityMapping.hpp"
#include "Minimizations/Mappings/LpDualityMapping.hpp"
#include "Minimizations/Mappings/LInfinityDualityMapping.hpp"

Mapping_ptr_t
PowerTypeDualityMappingFactory::createInstance(
		const NormedSpace_ptr_t &_NormedSpaceRef,
		const double _p,
		const double _power)
{
	PowerTypeDualityMapping *mapping = NULL;
	if (_p == std::numeric_limits<double>::infinity()) {
		mapping = new LInfinityDualityMapping(_NormedSpaceRef,_power);
	} else if (_p > 1.) {
		mapping = new LpDualityMapping(_NormedSpaceRef,_p,_power);
	} else if (_p == 1.) {
		mapping = new L1DualityMapping(_NormedSpaceRef,_power);
	}
	return Mapping_ptr_t(mapping);
}
