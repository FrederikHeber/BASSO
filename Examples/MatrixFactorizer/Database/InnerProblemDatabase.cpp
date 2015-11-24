/*
 * InnerProblemDatabase.cpp
 *
 *  Created on: Nov 23, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InnerProblemDatabase.hpp"

#include <algorithm>

#include <boost/assign.hpp>

#include "Database/Table.hpp"
#include "Log/Logging.hpp"
#include "MatrixFactorizer/Database/MaxMinAverage.hpp"

using namespace boost::assign;

InnerProblemDatabase::InnerProblemDatabase()
{}

InnerProblemDatabase::~InnerProblemDatabase()
{
	BOOST_LOG_TRIVIAL(info)
			<< "The directory contains " << size() << " tables.";
	for (TableDirectory::TableIterator_t iter = directory.begin();
			iter != directory.end(); ++iter)
		BOOST_LOG_TRIVIAL(info)
				<< "Table " << iter->first << " has " << iter->second->size() << " tuples.";
}

std::vector<Table::any_values_t> convertTuplesToValueVector(
		const Table &_table,
		const Table::KeyType_t &_KeyTypes,
		const Table::keys_t &_allowed_keys)
{
	std::vector<Table::any_values_t> valuevector;
	valuevector.reserve(_allowed_keys.size());

	Table::keys_t tempkeys(_allowed_keys);
	for (Table::KeyType_t::const_iterator keytypeiter = _KeyTypes.begin();
			keytypeiter != _KeyTypes.end(); ++keytypeiter) {
		if (tempkeys.find(keytypeiter->first) != tempkeys.end()) {
			Table::any_values_t values =
					_table.getAllValuesPerType(keytypeiter->first);
			valuevector.push_back(values);
			tempkeys.erase(keytypeiter->first);
		}
	}

	return valuevector;
}

template <class T>
std::vector<T> convertAnyVectorTo(const Table::any_values_t &_valuevector)
{
	std::vector<T> returnvector;
	returnvector.reserve(_valuevector.size());
	for (Table::any_values_t::const_iterator iter = _valuevector.begin();
			iter != _valuevector.end(); ++iter)
		returnvector.push_back(boost::get<T>(
				boost::any_cast<Database_types::typevariant_t>(*iter)));

	return returnvector;
}

void InnerProblemDatabase::accumulateData() const
{
	// clear all present data
	accumulatedData.clear();

	// fetch table
	const Table& datatable = getTableConst("data_per_iteration");

	// go through tables and fetch values
	const Table::keys_t keys = datatable.getSetofUniqueKeys();
	Table::KeyType_t KeyTypes = datatable.getKeyToTypeMap(keys);
	assert(keys.size() == KeyTypes.size());
	Table::keys_t allowed_keys;
	std::string keyname("residual");
	allowed_keys += keyname;
	assert(std::includes(keys.begin(), keys.end(), allowed_keys.begin(), allowed_keys.end()));
	std::vector<Table::any_values_t> valuevectors =
			convertTuplesToValueVector(datatable, KeyTypes, allowed_keys);
	assert(valuevectors.size() == 1);

	// accumulate values and store
	BOOST_LOG_TRIVIAL(info)
			<< "There are " << valuevectors[0].size() << " values for " << keyname;
	MaxMinAverage<double> Bregman(convertAnyVectorTo<double>(valuevectors[0]));
	accumulatedData["min_"+keyname] = Bregman.minimum;
	accumulatedData["max_"+keyname] = Bregman.maximum;
	accumulatedData["avg_"+keyname] = Bregman.average;
	accumulatedData["var_"+keyname] = Bregman.variance;
}

void InnerProblemDatabase::prepareTableForAccumulatedValues(
		Table &_table,
		const std::string &_keyname)
{
	Table::Tuple_t &tuple = _table.getTuple();

	DefaultValue<double> defaultvalue;
	tuple.insert( std::make_pair("min_"+_keyname, defaultvalue.get()), Table::Data);
	tuple.insert( std::make_pair("max_"+_keyname, defaultvalue.get()), Table::Data);
	tuple.insert( std::make_pair("avg_"+_keyname, defaultvalue.get()), Table::Data);
	tuple.insert( std::make_pair("var_"+_keyname, defaultvalue.get()), Table::Data);
}

void InnerProblemDatabase::insertAccumulatedValues(
		Table &_table) const
{
	accumulateData();

	Table::Tuple_t &tuple = _table.getTuple();

	for (Table::TokenTypeMap_t::const_iterator iter = accumulatedData.begin();
			iter != accumulatedData.end();
			++iter) {
		tuple.replace( iter->first, iter->second);
	}
}
