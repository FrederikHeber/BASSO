/*
 * NormedSpaceFactory.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormedSpaceFactory.hpp"

#include "Math/Helpers.hpp"
//#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"


NormedSpace_ptr_t NormedSpaceFactory::createInstance(
		const unsigned int _dimension,
		const double _p)
{
	// create two empty spaces
	NormedSpace_ptr_t instance(
			new NormedSpace(_dimension) );
	instance->setSpace( instance );
	NormedSpace_ptr_t dualinstance(
			new NormedSpace(_dimension) );
	dualinstance->setSpace( dualinstance );

	// and link the dual spaces
	instance->setDualSpace(dualinstance);
	dualinstance->setDualSpace(instance);

	// calculate conjugate value
	const double q = Helpers::ConjugateValue(_p);

	// create norm instances and hand over to spaces
	Norm_ptr_t norm = NormFactory::createLpInstance(_p);
	instance->setNorm(norm);
	Norm_ptr_t dualnorm = NormFactory::createLpInstance(q);
	dualinstance->setNorm(dualnorm);

	// create duality mapping instance
//	Mapping_ptr_t mapping(
//			new Mapping(instance, dualinstance, p));
//	instance->setDualityMapping(mapping);
//	Mapping_ptr_t dualmapping =
//			mapping->getAdjointMapping();
//	dualinstance->setDualityMapping(dualmapping);

	return instance;
}


