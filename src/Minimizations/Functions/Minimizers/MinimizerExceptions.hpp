/*
 * MinimizerExceptions.hpp
 *
 *  Created on: Jan 20, 2015
 *      Author: heber
 */

#ifndef MINIMIZEREXCEPTIONS_HPP_
#define MINIMIZEREXCEPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/exception/all.hpp>
#include <string>

//!> typedef for error information, here a string with the variable name
typedef boost::error_info<
		struct tag_MinimizerIllegalNumber_variablename,
		std::string> MinimizerIllegalNumber_variablename;

/** Exception structure when minimum argument is illegal (inf, nan).
 *
 */
struct MinimizerIllegalNumber_exception :
	virtual boost::exception,
	virtual std::exception
{ };

#endif /* MINIMIZEREXCEPTIONS_HPP_ */
