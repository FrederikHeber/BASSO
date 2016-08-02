/*
 * DatabaseExceptions.hpp
 *
 *  Created on: Aug 2, 2016
 *      Author: heber
 */

#ifndef DATABASEEXCEPTIONS_HPP_
#define DATABASEEXCEPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/exception/all.hpp>
#include <string>

//!> typedef for error information, here a string with the variable name
typedef boost::error_info<
		struct tag_DatabaseIllegalKey_name,
		std::string> DatabaseIllegalKey_name;

/** Exception structure for illegal values handed to a routine.
 *
 */
struct DatabaseIllegalKey_exception :
	virtual boost::exception,
	virtual std::exception
{ };

#endif /* DATABASEEXCEPTIONS_HPP_ */
