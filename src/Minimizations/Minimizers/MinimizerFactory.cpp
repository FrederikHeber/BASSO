/*
 * MinimizerFactory.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "MinimizerFactory.hpp"

#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/LandweberMinimizer.hpp"
#include "Minimizations/Minimizers/SequentialSubspaceMinimizer.hpp"
#include "Minimizations/Minimizers/SequentialSubspaceMinimizerNoise.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion.hpp"
#include "Options/CommandLineOptions.hpp"

// static entities
const std::string MinimizerFactory::TypeNames[] = {
		"Landweber",
		"SESOP",
		"RESESOP"
};

MinimizerFactory::instance_ptr_t
MinimizerFactory::createInstance(
		const CommandLineOptions &_opts,
		const InverseProblem_ptr_t &_inverseproblem,
		Database &_database,
		const StoppingCriterion::ptr_t &_stopping_criteria
		)
{
	// create the instance depending on the type
	GeneralMinimizer *instance = NULL;
	switch(_opts.type) {
	case landweber:
		instance = new LandweberMinimizer(
				_inverseproblem,
				_opts.delta,
				_opts.maxiter,
				_opts.maxinneriter,
				_database,
				_stopping_criteria,
				(const enum DetermineStepWidthFactory::stepwidth_enumeration)_opts.stepwidth_type,
				_opts.outputsteps
				);
		break;
	case sequentialsubspace:
			instance = new SequentialSubspaceMinimizer(
					_inverseproblem,
					_opts.delta,
					_opts.maxiter,
					_opts.maxinneriter,
					_database,
					_stopping_criteria,
					_opts.outputsteps,
					_opts.orthogonalization_type
					);
			break;
	case sequentialsubspace_noise:
			instance = new SequentialSubspaceMinimizerNoise(
					_inverseproblem,
					_opts.delta,
					_opts.maxiter,
					_opts.maxinneriter,
					_database,
					_stopping_criteria,
					_opts.outputsteps,
					_opts.orthogonalization_type
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
