/*
 * FunctionalMinimizerFactory.cpp
 *
 *  Created on: Dec 11, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "FunctionalMinimizerFactory.hpp"

#include <algorithm>
#include <iterator>

#include "Log/Logging.hpp"

// static entities
const std::string FunctionalMinimizerFactory::TypeNames[] = {
		"gsl",
		"nlopt"
};
enum FunctionalMinimizerFactory::MinimizationLibraries
FunctionalMinimizerFactory::CurrentMinLib(gnuscientificlibrary);
unsigned int FunctionalMinimizerFactory::maxiterations = 100;

bool FunctionalMinimizerFactory::isValidName(const std::string &_name)
{
	const std::string *begin = TypeNames;
	const std::string *end = TypeNames + MAX_MinimizationLibraries;
	return std::find(begin, end, _name) != &TypeNames[MAX_MinimizationLibraries];
}

void FunctionalMinimizerFactory::setMinLib(const std::string &_name)
{
	assert( isValidName(_name) );
	const std::string *begin = TypeNames;
	const std::string *end = TypeNames + MAX_MinimizationLibraries;
	CurrentMinLib = (enum MinimizationLibraries)std::distance(
			TypeNames,
			std::find( begin, end, _name));
}


