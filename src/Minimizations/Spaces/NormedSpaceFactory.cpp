/*
 * NormedSpaceFactory.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormedSpaceFactory.hpp"

#include <boost/bind.hpp>

#include "Math/Helpers.hpp"
#include "Minimizations/Mappings/IllegalDualityMapping.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMappingFactory.hpp"
#include "Minimizations/Mappings/Specifics/RelativeShrinkageMapping.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Norms/Specifics/RegularizedL1Norm.hpp"
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
			new NormedDualSpace(_dimension) );
	dualinstance->setSpace( dualinstance );

	// and link the (now two) dual spaces
	instance->setDualSpace(dualinstance);
	dualinstance->setDualSpace(instance);

	// create norm
	Norm_ptr_t norm = NormFactory::getInstance().create(
			"lp",
			instance,
			NormFactory::args_t(1, boost::any(_p)));
	instance->setNorm(norm);

	// create dual norm
	Norm_ptr_t dualnorm = NormFactory::getInstance().create(
			"dual_lp",
			dualinstance,
			NormFactory::args_t(1, boost::any(_p)));
	dualinstance->setNorm(dualnorm);

	// create duality mapping instance function call
	NormedSpace::constructDualityMapping_t mapping_cstor =
			boost::bind(
					&PowerTypeDualityMappingFactory::createInstance,
					NormedSpace_weakptr_t(instance),
					_power);
	instance->setDualityMappingConstructor(mapping_cstor);

	// create duality mapping instance constructor. We need this as both
	// mappings need to be created at the same time (and inter-linked,
	// similarly as with the spaces)
	// TODO: We could use a DualityMappingFactory for this instead: Then it
	// would be consistent with how spaces are constructed
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
			new NormedDualSpace(_dimension) );
	dualinstance->setSpace( dualinstance );

	// and link the (now two) dual spaces
	instance->setDualSpace(dualinstance);
	dualinstance->setDualSpace(instance);

	// create norm instances and hand over to spaces
	Norm_ptr_t norm = NormFactory::getInstance().create(
			"regularized_l1",
			instance,
			NormFactory::args_t(1, boost::any(_lambda)));
	instance->setNorm(norm);

	// create dual norm
	Norm_ptr_t dualnorm = NormFactory::getInstance().create(
			"dual_regularized_l1",
			instance,
			NormFactory::args_t(1, boost::any(_lambda)));
	dualinstance->setNorm(dualnorm);

	// create duality mapping instance: we only have the mapping from the
	// dual space back into the source space, not the other way round, as
	// the source space is not smooth (hence no single-valued duality
	// mapping exists).
	Mapping_ptr_t mapping( new IllegalDualityMapping );
	instance->setDualityMapping(mapping);

	// create duality mapping instance
	Mapping_ptr_t dualmapping(
			new RelativeShrinkageMapping(
					dualinstance,
					_lambda)
	);
	dualinstance->setDualityMapping(dualmapping);

	return instance;
}
