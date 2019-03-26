/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
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
		Database &_database
		)
{
	// create the instance depending on the type
	GeneralMinimizer *instance = NULL;
	switch(_opts.type) {
	case landweber:
		instance = new LandweberMinimizer(
				_opts,
				_inverseproblem,
				_database
				);
		break;
	case sequentialsubspace:
			instance = new SequentialSubspaceMinimizer(
					_opts,
					_inverseproblem,
					_database
					);
			break;
	case sequentialsubspace_noise:
			instance = new SequentialSubspaceMinimizerNoise(
					_opts,
					_inverseproblem,
					_database
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
