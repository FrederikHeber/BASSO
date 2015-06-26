/*
 * Logging.cpp
 *
 *  Created on: Jun 11, 2014
 *      Author: heber
 */

#include "Log/Logging.hpp"

#include "BassoConfig.h"

BOOST_LOG_ATTRIBUTE_KEYWORD(timestamp, "TimeStamp", boost::posix_time::ptime);

void startLogging()
{
	// change log format
	logging::add_console_log(std::cout,
			keywords::format = (
				expr::stream
					<< expr::format_date_time(timestamp, "%Y-%m-%d %H:%M:%S")
					//<< expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
					<< ": <" << logging::trivial::severity
					<< "> " << expr::smessage
					)
			//keywords::format = "[%TimeStamp%]: [%Severity%] %Message%"
			);
	// register attributes
	logging::add_common_attributes();
}

void stopLogging()
{
    boost::shared_ptr< logging::core > core = logging::core::get();

    // Remove the sink from the core, so that no records are passed to it
    core->remove_all_sinks();

    // flush output
    core->flush();
}

