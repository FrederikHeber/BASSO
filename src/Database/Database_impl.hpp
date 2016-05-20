/*
 * Database_impl.hpp
 *
 *  Created on: Dec 12, 2014
 *      Author: heber
 */

#ifndef DATABASE_DATABASE_IMPL_HPP_
#define DATABASE_DATABASE_IMPL_HPP_

#include "BassoConfig.h"

#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/repetition/enum_shifted_binary_params.hpp>
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <iterator>
#include <sstream>

#define QuestionmarkPrinter(z, count, data) \
	BOOST_PP_COMMA_IF(count) ?

#define QuestionmarksPrinter(count) \
	VALUES( BOOST_PP_REPEAT( count, QuestionmarkPrinter, ~) )

#define ComponentPrinter(z, n, vectorname) \
	, use(vectorname [ n ])

inline
std::string printVector(const std::vector<std::string> &_vector)
{
	std::stringstream stream;
	if (_vector.empty())
		return std::string("");
	else if (_vector.size() == 1)
		stream << _vector[0];
	else
		std::copy(_vector.begin(), _vector.end(),
				(std::ostream_iterator<std::string>(stream, ",")));
	return stream.str();
}

#define StringizedComponentPrinter(z, n, vectorname) \
	<< ", use( " << \
			printVector(vectorname[n]) \
	<< " )"

#define ComponentsPrinter(count, vectorname) \
	BOOST_PP_REPEAT( count, ComponentPrinter, vectorname)

#define StringizedComponentsPrinter(count, vectorname) \
	BOOST_PP_REPEAT( count, StringizedComponentPrinter, vectorname)

#define CasePrinter(z, count, argumentsequence) \
	case count : \
	{ \
		LOG(trace, "SQL: " \
		<< "INSERT INTO " << \
		BOOST_PP_SEQ_ELEM(1, argumentsequence) \
		<< " " << \
		BOOST_PP_STRINGIZE( QuestionmarksPrinter(count) ) \
		StringizedComponentsPrinter(count, \
				BOOST_PP_SEQ_ELEM(2, argumentsequence) \
				)); \
		BOOST_PP_SEQ_ELEM(0, argumentsequence) \
		<< "INSERT INTO " << \
		BOOST_PP_SEQ_ELEM(1, argumentsequence) \
		<< " " << \
		BOOST_PP_STRINGIZE( QuestionmarksPrinter(count) ) \
		ComponentsPrinter(count, \
				BOOST_PP_SEQ_ELEM(2, argumentsequence) \
				) \
		, now; \
		break; \
	}


#endif /* DATABASE_DATABASE_IMPL_HPP_ */
