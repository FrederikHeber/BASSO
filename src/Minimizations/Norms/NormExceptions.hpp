/*
 * NormExceptions.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef NORMEXCEPTIONS_HPP_
#define NORMEXCEPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/exception/all.hpp>
#include <string>

//!> typedef for error information, here a string with the variable name
typedef boost::error_info<
		struct tag_NormIllegalValue_name,
		std::string> NormIllegalValue_name;

/** Exception structure for illegal values handed to a routine.
 *
 */
struct NormIllegalValue_exception :
	virtual boost::exception,
	virtual std::exception
{ };

#endif /* NORMEXCEPTIONS_HPP_ */
