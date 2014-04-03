/*
 * MathExceptions.hpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#ifndef MATHEXCEPTIONS_HPP_
#define MATHEXCEPTIONS_HPP_

#include <boost/exception/all.hpp>
#include <string>

//!> typedef for error information, here a string with the variable name
typedef boost::error_info<
		struct tag_MathIllegalValue_name,
		std::string> MathIllegalValue_name;

/** Exception structure for illegal values handed to a routine.
 *
 */
struct MathIllegalValue_Error:
	virtual boost::exception,
	virtual std::exception
{ };

#endif /* MATHEXCEPTIONS_HPP_ */
