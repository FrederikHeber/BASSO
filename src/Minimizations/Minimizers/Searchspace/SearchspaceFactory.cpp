/*
 * SearchspaceFactory.cpp
 *
 *  Created on: Nov 18, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SearchspaceFactory.hpp"

#include <cassert>

// static entities
SearchspaceFactory::NameTypeMap_t SearchspaceFactory::NameTypeMap;

const std::string
SearchspaceFactory::getName(const enum SearchspaceType _type)
{
	fillNameTypeMap();
	if (isValidType(_type))
		return NameTypeMap.left.at(_type);
	else
		return NameTypeMap.left.at(InvalidType);
}

const enum SearchspaceFactory::SearchspaceType
SearchspaceFactory::getType(const std::string &_name)
{
	fillNameTypeMap();
	if (isValidName(_name))
		return NameTypeMap.right.at(_name);
	else
		return InvalidType;
}

bool
SearchspaceFactory::isValidName(const std::string &_name)
{
	fillNameTypeMap();
	NameTypeMap_t::right_const_iterator iter =
			NameTypeMap.right.find(_name);
	return ((iter != NameTypeMap.right.end())
			&& (iter->second != InvalidType));
}

bool
SearchspaceFactory::isValidType(const enum SearchspaceType _type)
{
	fillNameTypeMap();
	NameTypeMap_t::left_const_iterator iter =
			NameTypeMap.left.find(_type);
	return ((iter != NameTypeMap.left.end())
			&& (iter->first != InvalidType));
}

void
SearchspaceFactory::fillNameTypeMap()
{
	if (NameTypeMap.empty()) {
		NameTypeMap.insert(
				NameTypeMap_t::value_type( InvalidType, "") );
		NameTypeMap.insert(
				NameTypeMap_t::value_type( LastNDirections, "LastNDirections") );
		NameTypeMap.insert(
				NameTypeMap_t::value_type( Nemirovsky, "Nemirovsky") );
		// check that we inserted all
		assert( NameTypeMap.size() == (size_t)MAX_SearchspaceType );
	}
}
