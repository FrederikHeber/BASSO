/*
 * MinimizerFactory.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "MinimizerFactory.hpp"

#include "Minimizations/DualityMappings/DualityMappingsContainer.hpp"
#include "Minimizations/DualityMappings/DefaultDualityMappings.hpp"
#include "Minimizations/DualityMappings/RegularizedL1Norm.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/LandweberMinimizer.hpp"
#include "Minimizations/Minimizers/SequentialSubspaceMinimizer.hpp"
#include "Minimizations/Minimizers/SequentialSubspaceMinimizerNoise.hpp"

// static entities
const std::string MinimizerFactory::TypeNames[] = {
		"landweber",
		"SSO",
		"SSO_noise"
};

DualityMappingsContainer *MinimizerFactory::DualityContainer = NULL;

MinimizerFactory::instance_ptr_t MinimizerFactory::createInstance(
		const enum InstanceType &_type,
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		Database &_database,
		const unsigned int _outputsteps
		)
{
	// create DualityMappingsContainer instance
	if (DualityContainer != NULL)
		delete DualityContainer;
	DualityContainer = new DefaultDualityMappings(
			_inverseproblem->x->getSpace()->getNorm()->getPvalue(),
			_inverseproblem->x->getSpace()->getDualityMapping()->getPower(),
			1e-6
			);

	// return the instance depending on the type
	return getMinimizerInstance(
					_type,
					_inverseproblem,
					_Delta,
					_maxiter,
					_database,
					_outputsteps
			);
}

MinimizerFactory::instance_ptr_t MinimizerFactory::getRegularizedInstance(
		const enum InstanceType &_type,
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		Database &_database,
		const unsigned int _outputsteps
		)
{
	// create DualityMappingsContainer instance
	if (DualityContainer != NULL)
		delete DualityContainer;
	const Mapping_ptr_t &J_q = _inverseproblem->y->getSpace()->getDualityMapping();
	const SoftThresholdingMapping & S =
			dynamic_cast<SoftThresholdingMapping &>(*J_q);
	DualityContainer = new RegularizedL1Norm(
			S.getLambda()
			);

	// return the instance depending on the type
	return getMinimizerInstance(
					_type,
					_inverseproblem,
					_Delta,
					_maxiter,
					_database,
					_outputsteps
			);
}


MinimizerFactory::instance_ptr_t
MinimizerFactory::getMinimizerInstance(
		const enum InstanceType &_type,
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		Database &_database,
		const unsigned int _outputsteps
		)
{
	// create the instance depending on the type
	GeneralMinimizer *instance = NULL;
	switch(_type) {
	case landweber:
		instance = new LandweberMinimizer(
				_inverseproblem,
				_Delta,
				_maxiter,
				_database,
				_outputsteps
				);
		break;
	case sequentialsubspace:
			instance = new SequentialSubspaceMinimizer(
					_inverseproblem,
					_Delta,
					_maxiter,
					_database,
					_outputsteps
					);
			break;
	case sequentialsubspace_noise:
			instance = new SequentialSubspaceMinimizerNoise(
					_inverseproblem,
					_Delta,
					_maxiter,
					_database,
					_outputsteps
					);
			break;
	default:
		std::cerr << "Illegal or unknown type of GeneralMinimizer requested."
			<< std::endl;
		break;
	}

	// return the wrapped instance
	return instance_ptr_t(instance);
}



unsigned int MinimizerFactory::getTypeNamesIndex(
		const std::string &_name)
{
	unsigned int i=0;
	for (; i<MAX_InstanceType; ++i)
		if (TypeNames[i] == _name)
			break;

	return i;
}

MinimizerFactory::InstanceType MinimizerFactory::getTypeForName(
		const std::string &_name)
{
	unsigned int i=getTypeNamesIndex(_name);
	assert(i != MAX_InstanceType);
	return (enum InstanceType)i;
}

const std::string& MinimizerFactory::getNameForType(
		const enum InstanceType &_type)
{
	assert( (_type >= landweber) && (_type < MAX_InstanceType) );
	return TypeNames[_type];
}
