/*
 * NormedSpaceFactory.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormedSpaceFactory.hpp"

#include "Math/Helpers.hpp"
#include "Minimizations/Mappings/IllegalDualityMapping.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMappingFactory.hpp"
#include "Minimizations/Mappings/SoftThresholdingMapping.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

// static instance
const NormedSpace_ptr_t NormedSpaceFactory::DummySpace =
		NormedSpaceFactory::createLpInstance(
				0, 2., 2.);

NormedSpace_ptr_t NormedSpaceFactory::createLpInstance(
		const unsigned int _dimension,
		const double _p,
		const double _power)
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
	Norm_ptr_t norm = NormFactory::createLpInstance(instance, _p);
	instance->setNorm(norm);
	Norm_ptr_t dualnorm = NormFactory::createLpInstance(dualinstance, q);
	dualinstance->setNorm(dualnorm);

	// create duality mapping instance
	Mapping_ptr_t mapping =
			PowerTypeDualityMappingFactory::createInstance(
					instance, _power);
	instance->setDualityMapping(mapping);
	Mapping_ptr_t dualmapping =
			mapping->getAdjointMapping();
	dualinstance->setDualityMapping(dualmapping);

	return instance;
}

NormedSpace_ptr_t NormedSpaceFactory::createRegularizedL1Instance(
		const unsigned int _dimension,
		const double _lambda,
		const double _power)
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
//	const double q = std::numeric_limits<double>::infinity();

	// create norm instances and hand over to spaces
	Norm_ptr_t norm = NormFactory::createRegularizedL1Instance(
			instance, _lambda);
	instance->setNorm(norm);
	// TODO: I don't know the dual norm to regularized l1 norm yet
	Norm_ptr_t dualnorm = NormFactory::createIllegalInstance(dualinstance);
	dualinstance->setNorm(dualnorm);

	// create duality mapping instance: we only have the mapping from
	// the dual space into the source space, not the other way round as
	// the source space is not smooth (hence not single-valued duality
	// mapping exists).
	Mapping_ptr_t mapping(
				new IllegalDualityMapping
	);
	instance->setDualityMapping(mapping);
	Mapping_ptr_t dualmapping(
			new SoftThresholdingMapping(
					dualinstance,
					_lambda)
	);
	dualinstance->setDualityMapping(dualmapping);

	return instance;
}
