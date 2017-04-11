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
 * Logging.cpp
 *
 *  Created on: Jun 11, 2014
 *      Author: heber
 */

#include "Log/Logging.hpp"

#include "BassoConfig.h"

#include "version.h"

BOOST_LOG_ATTRIBUTE_KEYWORD(timestamp, "TimeStamp", boost::posix_time::ptime);

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace keywords = boost::log::keywords;


void startLogging()
{
	// change log format
	boost::log::add_console_log(std::cout,
			keywords::format = (
				expr::stream
					<< expr::format_date_time(timestamp, "%Y-%m-%d %H:%M:%S")
					//<< expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
					<< ": <" << boost::log::trivial::severity
					<< "> " << expr::smessage
					)
			//keywords::format = "[%TimeStamp%]: [%Severity%] %Message%"
			);
	// register attributes
	boost::log::add_common_attributes();
}

void stopLogging()
{
    boost::shared_ptr< boost::log::core > core = boost::log::core::get();

    // Remove the sink from the core, so that no records are passed to it
    core->remove_all_sinks();

    // flush output
    core->flush();
}

void showVersion(const std::string _programname)
{
	if (_programname.empty())
		std::cout << "Version " << Basso_VERSION_MAJOR << "." << Basso_VERSION_MINOR
			<< " -- build " << VERSION << std::endl;
	else
		std::cout << _programname << " version " << Basso_VERSION_MAJOR << "."
		<< Basso_VERSION_MINOR << " -- build " << VERSION << std::endl;
}

void showCopyright()
{
	std::cout << "(C) 2014-2016 Universität des Saarlandes. All rights reserved."
			<< std::endl;
}

