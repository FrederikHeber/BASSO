/*
 * MinimizerFactory.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "Minimizations/GeneralMinimizer.hpp"
#include "Minimizations/LandweberMinimizer.hpp"
#include "Minimizations/MinimizerFactory.hpp"
#include "Minimizations/SequentialSubspaceMinimizer.hpp"


// static entities
const std::string MinimizerFactory::TypeNames[] = {
		"landweber",
		"SSO"
};

MinimizerFactory::instance_ptr_t
MinimizerFactory::getInstance(
		const enum InstanceType &_type,
		const double _NormX,
		const double _NormY,
		const double _PowerX,
		const double _PowerY,
		const double _Delta,
		const unsigned int _maxiter,
		const unsigned int _outputsteps
		)
{
	// create the instance depending on the type
	GeneralMinimizer *instance = NULL;
	switch(_type) {
	case landweber:
		instance = new LandweberMinimizer(
				_NormX,
				_NormY,
				_PowerX,
				_PowerY,
				_Delta,
				_maxiter,
				_outputsteps
				);
		break;
	case sequentialsubspace:
			instance = new SequentialSubspaceMinimizer(
					_NormX,
					_NormY,
					_PowerX,
					_PowerY,
					_Delta,
					_maxiter,
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
