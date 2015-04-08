/*
 * types.hpp
 *
 *  Created on: Apr 8, 2015
 *      Author: heber
 */

#ifndef DATABASE_TYPES_HPP_
#define DATABASE_TYPES_HPP_

#include "BassoConfig.h"

#include <boost/variant.hpp>

namespace Database_types {
	//!> defines a variant of all types essential to a database
	typedef boost::variant<int, double, std::string > typevariant_t;

	//!> enumeration of all type variants used as database types
	enum types_t
	{
		inttype=0,
		doubletype=1,
		valchartype=2,
		MAX_TYPES
	};
} /* End namespace Database_types */

#endif /* DATABASE_TYPES_HPP_ */
