/*
 * Table_Tuple.cpp
 *
 *  Created on: Aug 22, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "Table.hpp"

#include <cassert>
#include <string>

#include "Database/types.hpp"

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
		const Database_types::typevariant_t &_value)
{
	std::pair<iterator, bool> inserter =
			TokenTypeMap_t::insert(
							std::make_pair(_key, _value) );
	assert( !inserter.second );
	inserter.first->second = _value;
}

void Table::Tuple_t::insert(const std::pair<std::string,
		Database_types::typevariant_t> &_pair,
		const enum ColumnType _type)
{
	{
		std::pair<iterator, bool> inserter =
				TokenTypeMap_t::insert(
								_pair );
		assert( inserter.second );
	}
	{
		std::pair<TypeMap_t::iterator, bool> inserter =
				TypeMap.insert( std::make_pair(_pair.first, _type) );
		assert( inserter.second );
	}
}

bool Table::Tuple_t::isParameter(const std::string &_name) const
{
	TypeMap_t::const_iterator iter = TypeMap.find(_name);
	assert( iter != TypeMap.end() );
	return (iter->second == Parameter);
}

void Table::Tuple_t::clear()
{
	TokenTypeMap_t::clear();
	TypeMap.clear();
}
