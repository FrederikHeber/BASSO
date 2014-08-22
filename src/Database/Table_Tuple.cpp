/*
 * Table_Tuple.cpp
 *
 *  Created on: Aug 22, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "Table.hpp"

#include <string>

#include "Database.hpp"

bool Table::Tuple_t::operator<(const Table::Tuple_t &_a) const
{
	bool status = true;
	Tuple_t::const_iterator firstiter = begin();
	Tuple_t::const_iterator seconditer = _a.begin();
	for (;(firstiter != end()) && (seconditer != _a.end());
		++firstiter, ++seconditer) {
		if (firstiter->first < seconditer->first)
			break;
		else if (firstiter->first > seconditer->first) {
			status = false;
			break;
		} else {
			// check values
		}
	}
	return status;
}

void Table::Tuple_t::replace(
		const std::string &_key,
		const Database::typevariant_t &_value)
{
	std::pair<iterator, bool> inserter =
			insert( std::make_pair(_key, _value) );
	if (!inserter.second) {
		inserter.first->second = _value;
	}
}
