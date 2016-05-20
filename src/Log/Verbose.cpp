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

