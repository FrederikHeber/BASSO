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
#include "Minimizations/Mappings/Specifics/SoftThresholdingMapping.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Norms/RegularizedL1Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedDualSpace.hpp"

// static instance
const NormedSpace_ptr_t NormedSpaceFactory::DummySpace =
		NormedSpaceFactory::createLpInstance(
				0, 2., 2.);

NormedSpace_ptr_t NormedSpaceFactory::createLpInstance(
		const unsigned int _dimension,
		const double _p,
		const double _power)
{
	// create empty space
	NormedSpace_ptr_t instance(
			new NormedSpace(_dimension) );
	instance->setSpace( instance );
	NormedSpace_ptr_t dualinstance(
			new NormedDualSpace(instance->getDimension()) );
	dualinstance->setSpace( dualinstance );

	// and link the (now two) dual spaces
	instance->setDualSpace(dualinstance);
	dualinstance->setDualSpace(instance);

	// create norm
	Norm_ptr_t norm = NormFactory::createLpInstance(instance, _p);
	instance->setNorm(norm);

	// create dual norm
	const double q =
			Helpers::ConjugateValue(instance->getNorm()->getPvalue());
	Norm_ptr_t dualnorm = NormFactory::createLpInstance(dualinstance, q);
	dualinstance->setNorm(dualnorm);

	// create duality mapping instance function call
	NormedSpace::constructDualityMapping_t mapping_cstor =
			boost::bind(
					&PowerTypeDualityMappingFactory::createInstance,
					NormedSpace_weakptr_t(instance),
					_power);
	instance->setDualityMappingConstructor(mapping_cstor);

	// create duality mapping instance constructor
	const double qpower =
			Helpers::ConjugateValue(instance->getDualityMapping()->getPower());
	NormedSpace::constructDualityMapping_t dualmapping_cstor =
			boost::bind(
					&PowerTypeDualityMappingFactory::createInstance,
					/* only weakptr instance because is stored in this
					 * bound function and prevents shared_ptr from
					 * deallocating.
					 */
					NormedSpace_weakptr_t(dualinstance),
					qpower);
	dualinstance->setDualityMappingConstructor(dualmapping_cstor);

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
			new NormedDualSpace(instance->getDimension()) );
	dualinstance->setSpace( dualinstance );

	// and link the (now two) dual spaces
	instance->setDualSpace(dualinstance);
	dualinstance->setDualSpace(instance);

	// create norm instances and hand over to spaces
	Norm_ptr_t norm = NormFactory::createRegularizedL1Instance(
			instance, _lambda);
	instance->setNorm(norm);

	RegularizedL1Norm *regularized_norm =
			dynamic_cast<RegularizedL1Norm *>(instance->getNorm().get());
	// create dual norm
	Norm_ptr_t dualnorm =
			NormFactory::createDualRegularizedL1Instance(
					dualinstance, regularized_norm->getLambda());
	dualinstance->setNorm(dualnorm);

	// create duality mapping instance: we only have the mapping from
	// the dual space into the source space, not the other way round as
	// the source space is not smooth (hence not single-valued duality
	// mapping exists).
	Mapping_ptr_t mapping(
				new IllegalDualityMapping
	);
	instance->setDualityMapping(mapping);

	// create duality mapping instance
	Mapping_ptr_t dualmapping(
			new SoftThresholdingMapping(
					dualinstance,
					regularized_norm->getLambda())
	);
	dualinstance->setDualityMapping(dualmapping);

	return instance;
}
