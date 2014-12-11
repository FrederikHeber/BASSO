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
	/** Factory function creating a power type duality mapping, i.e. for
	 *  a Lp norm.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _power power type of the weight function
	 * @return
	 */
	static const Mapping_ptr_t createInstance(
			const NormedSpace_ptr_t &_NormedSpaceRef,
			const double _power);
};



#endif /* POWERTYPEDUALITYMAPPINGFACTORY_HPP_ */
