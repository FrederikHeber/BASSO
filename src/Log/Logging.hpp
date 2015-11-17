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

#include <iosfwd>
#include <iterator>
#include <list>
#include <set>
#include <vector>

void startLogging();
void stopLogging();

void showVersion(const std::string _programname = std::string(""));
void showCopyright();

template<class T>
std::ostream & operator<<(std::ostream &ost, const std::list<T> & values)
{
	std::copy(
			values.begin(), values.end(),
			std::ostream_iterator<T>(ost, " "));
	return ost;
}

template<class T>
std::ostream & operator<<(std::ostream &ost, const std::set<T> & values)
{
	std::copy(
			values.begin(), values.end(),
			std::ostream_iterator<T>(ost, " "));
	return ost;
}

template<class T>
std::ostream & operator<<(std::ostream &ost, const std::vector<T> & values)
{
	std::copy(
			values.begin(), values.end(),
			std::ostream_iterator<T>(ost, " "));
	return ost;
}

#endif /* LOGGING_HPP_ */
