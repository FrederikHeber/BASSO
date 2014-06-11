/*
 * Logging.hpp
 *
 *  Created on: Jun 11, 2014
 *      Author: heber
 */

#ifndef LOGGING_HPP_
#define LOGGING_HPP_

#include "BassoConfig.h"

#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace keywords = boost::log::keywords;

void startLogging();

#endif /* LOGGING_HPP_ */
