/*
 * DualityMappingFactory.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef DUALITYMAPPINGFACTORY_HPP_
#define DUALITYMAPPINGFACTORY_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

struct DualityMappingFactory
{
	/** Factory function creating a power type duality mapping, i.e. for
	 *  a Lp norm.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _power power type of the weight function
	 * @return
	 */
	static const Mapping_ptr_t createInstance(
			const NormedSpace_weakptr_t _NormedSpaceRef,
			const double _power);
};



#endif /* DUALITYMAPPINGFACTORY_HPP_ */
