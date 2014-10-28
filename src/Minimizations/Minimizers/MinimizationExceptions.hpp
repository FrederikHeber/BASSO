/*
 * MinimizationExceptions.hpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONEXCEPTIONS_HPP_
#define MINIMIZATIONEXCEPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/exception/all.hpp>
#include <string>

//!> typedef for error information, here a string with the variable name
typedef boost::error_info<
		struct tag_MinimizationIllegalValue_name,
		std::string> MinimizationIllegalValue_name;

//!> typedef for error information, here a function value (double)
typedef boost::error_info<
		struct tag_MinimizationFunctionError_name,
		double> MinimizationFunctionError_name;

/** Exception structure for illegal values handed to a routine.
 *
 */
struct MinimizationIllegalValue_exception :
	virtual boost::exception,
	virtual std::exception
{ };

/** Exception structure on function evaluation.
 *
 */
struct MinimizationFunctionError_exception :
	virtual boost::exception,
	virtual std::exception
{ };

#endif /* MINIMIZATIONEXCEPTIONS_HPP_ */
