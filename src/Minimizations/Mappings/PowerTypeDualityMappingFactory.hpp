/*
 * PowerTypeDualityMappingFactory.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef POWERTYPEDUALITYMAPPINGFACTORY_HPP_
#define POWERTYPEDUALITYMAPPINGFACTORY_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

struct PowerTypeDualityMappingFactory
{
	static PowerTypeDualityMapping_ptr_t createInstance(const double _p);
};



#endif /* POWERTYPEDUALITYMAPPINGFACTORY_HPP_ */
