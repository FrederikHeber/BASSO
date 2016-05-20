/*
 * Verbose.hpp
 *
 *  Created on: May 20, 2016
 *      Author: heber
 */

#ifndef LOG_VERBOSE_HPP_
#define LOG_VERBOSE_HPP_

#include "BassoConfig.h"

#include <cassert>
#include <vector>

#include <boost/log/trivial.hpp>

/** Contains the verbosity levels of the current execution.
 *
 */
struct Verbose
{
	/** Sets the desired level of verbosity.
	 *
	 * @param _level desired level
	 */
	static void setVerbosity(const unsigned int &_level);

	/** Getter for whether the specific \a _level of verbosity is
	 * admitted or not.
	 *
	 * @param _level specific level
	 * @return true - do output log, false - don't
	 */
	static bool verbosity(const enum boost::log::trivial::severity_level &_level) {
		assert( level_bool.size() == 6 );
		return level_bool[_level];
	}

	//!> typedef for level to bool mapping
	typedef std::vector<bool> level_bool_t;

private:

	//!> internal map containing the boolean result for each known level
	static level_bool_t level_bool;
};


#endif /* LOG_VERBOSE_HPP_ */
