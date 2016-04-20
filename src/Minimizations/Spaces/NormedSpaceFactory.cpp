/*
 * NormedSpaceFactory.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormedSpaceFactory.hpp"

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

#include "Math/Helpers.hpp"
#include <Minimizations/Mappings/DualityMappingFactory.hpp>
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedDualSpace.hpp"

using namespace boost::assign;

const NormedSpace_ptr_t NormedSpaceFactory::getDummySpace()
{
	static const NormedSpace_ptr_t DummySpace;
	if (DummySpace.get() == NULL) {
		NormedSpaceFactory::args_t args;
		args += boost::any(2.), boost::any(2.);
		const_cast<NormedSpace_ptr_t &>(DummySpace) =
				NormedSpaceFactory::create(
						0, "lp", args);
	}
	return DummySpace;
}

NormedSpace_ptr_t NormedSpaceFactory::create(
		const unsigned int _dimension,
		const std::string &_type,
		const args_t &_args)
{
	const std::string dualtype = "dual_"+_type;
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
	{
		Norm_ptr_t norm = NormFactory::getInstance().create(
				_type,
				instance,
				NormFactory::args_t(1, _args[0]));
		instance->setNorm(norm);
	}

	// create dual norm
	{
		Norm_ptr_t dualnorm = NormFactory::getInstance().create(
				dualtype,
				dualinstance,
				NormFactory::args_t(1, _args[0]));
		dualinstance->setNorm(dualnorm);
	}

	// create duality mapping instance
	{
		Mapping_ptr_t mapping =
				DualityMappingFactory::create(
						_type,
						NormedSpace_weakptr_t(instance),
						DualityMappingFactory::args_t(1, _args[1]));
		instance->setDualityMapping(mapping);
	}

	// create dual duality mapping instance
	{
		Mapping_ptr_t dualmapping =
				DualityMappingFactory::create(
						dualtype,
						NormedSpace_weakptr_t(dualinstance),
						DualityMappingFactory::args_t(1, _args[1]));
		dualinstance->setDualityMapping(dualmapping);
	}

	return instance;
}
