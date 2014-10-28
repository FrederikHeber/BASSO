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

PowerTypeDualityMapping_ptr_t
PowerTypeDualityMappingFactory::createInstance(
		const double _p,
		const double _power)
{
	PowerTypeDualityMapping *mapping = NULL;
	if (_p == std::numeric_limits<double>::infinity()) {
		mapping = new LInfinityDualityMapping(_power);
	} else if (_p > 1.) {
		mapping = new LpDualityMapping(_p,_power);
	} else if (_p == 1.) {
		mapping = new L1DualityMapping(_power);
	}
	return PowerTypeDualityMapping_ptr_t(mapping);
}