/*
 * SearchspaceFactory.cpp
 *
 *  Created on: Nov 18, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SearchspaceFactory.hpp"

#include <cassert>

#include "Log/Logging.hpp"

#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp"
#include "Minimizations/Minimizers/Searchspace/NemirovskyDirection.hpp"

// static entities
SearchspaceFactory::NameTypeMap_t SearchspaceFactory::NameTypeMap;
enum SearchspaceFactory::SearchspaceType
SearchspaceFactory::InstanceType = SearchspaceFactory::LastNDirections;

Searchspace::ptr_t SearchspaceFactory::create(
		const NormedSpace_ptr_t &_SearchDirectionSpace_ptr,
		const unsigned int _N,
		const LastNSearchDirections::OrthogonalizationType _orthogonalization_type
		)
#ifdef USE_OPENMP
#pragma omp critical
#endif /*USE_OPENMP */
{
	Searchspace::ptr_t returninstance;
	switch (InstanceType) {
	case LastNDirections:
		returninstance.reset(new LastNSearchDirections(
				_SearchDirectionSpace_ptr,
				_N,
				_orthogonalization_type)
				);
		break;
	case Nemirovsky:
		if (_N != 2) {
			LOG(error, "NemirovskyDirection always uses two search directions.");
			throw MinimizationIllegalValue_exception()
					<< MinimizationIllegalValue_name("N");
		}
		returninstance.reset(new NemirovskyDirection(
				_SearchDirectionSpace_ptr)
				);
		break;
	default:
		LOG(error, "Current InstanceType is unknown.");
		assert(0);
	}
	return returninstance;
}

void
SearchspaceFactory::setCurrentType(const enum SearchspaceType _type)
{
	if (isValidType(_type))
		InstanceType = _type;
}

const std::string
SearchspaceFactory::getName(const enum SearchspaceType _type)
#ifdef USE_OPENMP
#pragma omp critical
#endif /*USE_OPENMP */
{
	fillNameTypeMap();
	if (isValidType(_type))
		return NameTypeMap.left.at(_type);
	else
		return NameTypeMap.left.at(InvalidType);
}

const enum SearchspaceFactory::SearchspaceType
SearchspaceFactory::getType(const std::string &_name)
#ifdef USE_OPENMP
#pragma omp critical
#endif /*USE_OPENMP */
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
#ifdef USE_OPENMP
#pragma omp critical
#endif /*USE_OPENMP */
{
	fillNameTypeMap();
	NameTypeMap_t::left_const_iterator iter =
			NameTypeMap.left.find(_type);
	return ((iter != NameTypeMap.left.end())
			&& (iter->first != InvalidType));
}

void
SearchspaceFactory::fillNameTypeMap()
#ifdef USE_OPENMP
#pragma omp critical
#endif /*USE_OPENMP */
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
