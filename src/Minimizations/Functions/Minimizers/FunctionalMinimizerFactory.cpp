/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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


