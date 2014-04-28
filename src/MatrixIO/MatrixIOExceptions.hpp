/*
 * MatrixIOExceptions.hpp
 *
 *  Created on: Apr 28, 2014
 *      Author: heber
 */

#ifndef MATRIXIOEXCEPTIONS_HPP_
#define MATRIXIOEXCEPTIONS_HPP_

#include <boost/exception/all.hpp>
#include <string>

//!> typedef for error information, here a string with the file name
typedef boost::error_info<
		struct tag_MatrixIOError_filename,
		std::string> MatrixIOError_fileename;

/** Exception structure for stream ending too early.
 *
 */
struct MatrixIOStreamEnded_exception :
	virtual boost::exception,
	virtual std::exception
{ };

#endif /* MATRIXIOEXCEPTIONS_HPP_ */
