/*
 * PowerTypeDualityMappingFactory.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "PowerTypeDualityMappingFactory.hpp"

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Mappings/L1DualityMapping.hpp"
#include "Minimizations/Mappings/LpDualityMapping.hpp"
#include "Minimizations/Mappings/LInfinityDualityMapping.hpp"

PowerTypeDualityMapping_ptr_t
PowerTypeDualityMappingFactory::createInstance(const double _p)
{
	PowerTypeDualityMapping *mapping = NULL;
	if (_p > 1.) {
		mapping = new LpDualityMapping(_p);
	} else if (_p == 0.) {
		mapping = new LInfinityDualityMapping;
	} else if (_p == 1.) {
		mapping = new L1DualityMapping;
	}
	return PowerTypeDualityMapping_ptr_t(mapping);
}
