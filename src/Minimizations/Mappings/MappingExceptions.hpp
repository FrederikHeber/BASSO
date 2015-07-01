/*
 * MappingExceptions.hpp
 *
 *  Created on: Jul 01, 2015
 *      Author: heber
 */

#ifndef MAPPINGEXCEPTIONS_HPP_
#define MAPPINGEXCEPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/exception/all.hpp>
#include <string>

//!> typedef for error information, here a string with the variable name
typedef boost::error_info<
		struct tag_MappingIllegalValue_name,
		std::string> MappingIllegalValue_name;

/** Exception structure for illegal values handed to a routine.
 *
 */
struct MappingIllegalValue_exception :
	virtual boost::exception,
	virtual std::exception
{ };

#endif /* MAPPINGEXCEPTIONS_HPP_ */
