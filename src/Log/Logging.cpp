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

