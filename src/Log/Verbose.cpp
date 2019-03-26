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
 * Verbose.cpp
 *
 *  Created on: May 20, 2016
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Verbose.hpp"

#include "Logging.hpp"

static Verbose::level_bool_t initialVerbosity()
{
	Verbose::level_bool_t returnvalues(6, true);
	returnvalues[boost::log::trivial::trace] = false;
	returnvalues[boost::log::trivial::debug] = false;
	return returnvalues;
}

// static instances
Verbose::level_bool_t Verbose::level_bool(initialVerbosity());

void Verbose::setVerbosity(const unsigned int &_level)
{
	// stop (old) logging
	stopLogging();

	// set internal map
	level_bool.clear();
	level_bool.resize(6, true);

	// set filters in boost::log
	switch (_level) {
	default:
	case 0:
		for (unsigned int i=0;i<=1;++i)
			level_bool[i] = false;
		boost::log::core::get()->set_filter
		(
				boost::log::trivial::severity >= boost::log::trivial::info
		);
		break;
	case 1:
		level_bool[boost::log::trivial::trace] = false;
		boost::log::core::get()->set_filter
		(
				boost::log::trivial::severity >= boost::log::trivial::debug
		);
		break;
	case 2:
		boost::log::core::get()->set_filter
		(
				boost::log::trivial::severity >= boost::log::trivial::trace
		);
		break;
	}

	// start logging
	startLogging();
}

