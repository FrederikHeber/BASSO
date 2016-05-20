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

#include "Verbose.hpp"

#define LOG(level, msg) \
	{ if (Verbose::verbosity(boost::log::trivial:: level)) \
		BOOST_LOG_TRIVIAL(level) << msg; }

//!> typedef for the long formatting_ostream of boost::log
typedef boost::log::v2s_mt_posix::formatting_ostream boost_log_ostream;

void startLogging();
void stopLogging();

void showVersion(const std::string _programname = std::string(""));
void showCopyright();

template<class T>
inline
std::ostream & operator<<(std::ostream &ost, const std::list<T> & values)
{
	std::copy(
			values.begin(), values.end(),
			std::ostream_iterator<T>(ost, " "));
	return ost;
}

template<class T>
inline
std::ostream & operator<<(std::ostream &ost, const std::set<T> & values)
{
	std::copy(
			values.begin(), values.end(),
			std::ostream_iterator<T>(ost, " "));
	return ost;
}

template<class T>
inline
std::ostream & operator<<(std::ostream &ost, const std::vector<T> & values)
{
	std::copy(
			values.begin(), values.end(),
			std::ostream_iterator<T>(ost, " "));
	return ost;
}


template<class T>
inline
boost_log_ostream&
operator<<( boost_log_ostream&ost, const std::list<T> & values)
{
	for (typename std::list<T>::const_iterator iter = values.begin();
			iter != values.begin(); ++iter)
		ost << *iter << " ";
	return ost;
}

template<class T>
inline
boost_log_ostream&
operator<<( boost_log_ostream&ost, const std::set<T> & values)
{
	for (typename std::set<T>::const_iterator iter = values.begin();
			iter != values.begin(); ++iter)
		ost << *iter << " ";
	return ost;
}

template<class T>
inline
boost_log_ostream&
operator<<( boost_log_ostream&ost, const std::vector<T> & values)
{
	for (typename std::vector<T>::const_iterator iter = values.begin();
			iter != values.begin(); ++iter)
		ost << *iter << " ";
	return ost;
}

#endif /* LOGGING_HPP_ */
