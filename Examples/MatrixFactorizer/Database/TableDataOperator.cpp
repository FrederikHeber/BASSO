/*
 * TableDataOperator.cpp
 *
 *  Created on: Nov 24, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "TableDataOperator.hpp"

#include <algorithm>

#include "Database/Database.hpp"
#include "Log/Logging.hpp"

TableDataOperator::TableDataOperator(
		const Database&_db) :
	db(_db)
{}

TableDataOperator::TableDataOperator(
		const Database&_db,
		const std::string &_tablename,
		const keys_t::const_iterator &_accumulated_keys_begin,
		const keys_t::const_iterator &_accumulated_keys_end) :
	db(_db),
	tablename(_tablename)
{
	accumulated_keys.insert(_accumulated_keys_begin, _accumulated_keys_end);
}

TableDataOperator::valuevectors_t convertTuplesToValueVector(
		const Table &_table,
		const Table::KeyType_t &_KeyTypes,
		const Table::keys_t &_allowed_keys)
{
	TableDataOperator::valuevectors_t valuevector;
	valuevector.reserve(_allowed_keys.size());

	Table::keys_t tempkeys(_allowed_keys);
	for (Table::KeyType_t::const_iterator keytypeiter = _KeyTypes.begin();
			keytypeiter != _KeyTypes.end(); ++keytypeiter) {
		if (tempkeys.find(keytypeiter->first) != tempkeys.end()) {
			Table::any_values_t values =
					_table.getAllValuesPerType(keytypeiter->first);
			valuevector.push_back(
					std::make_pair(keytypeiter->second, values ));
			tempkeys.erase(keytypeiter->first);
		}
	}

	return valuevector;
}

TableDataOperator::valuevectors_t
TableDataOperator::extractValues() const
{
	// do nothing on empty keys
	if (accumulated_keys.empty())
		return valuevectors_t();

	// clear all present data
	accumulatedData.clear();

	// fetch table
	const Table& datatable = db.getTableConst(tablename);

	// go through tables and fetch values
	const Table::keys_t keys = datatable.getSetofUniqueKeys();
	{
		std::stringstream output;
		std::copy(keys.begin(), keys.end(),
				std::ostream_iterator<std::string>(output, ","));
		BOOST_LOG_TRIVIAL(trace)
			<< "keys are " << output.str();
	}
	{
		std::stringstream output;
		std::copy(accumulated_keys.begin(), accumulated_keys.end(),
				std::ostream_iterator<std::string>(output, ","));
		BOOST_LOG_TRIVIAL(trace)
			<< "accumulated_keys are " << output.str();
	}
	Table::KeyType_t KeyTypes = datatable.getKeyToTypeMap(keys);
	assert(keys.size() == KeyTypes.size());
	assert(std::includes(
			keys.begin(), keys.end(),
			accumulated_keys.begin(), accumulated_keys.end()));
	valuevectors_t valuevectors =
			convertTuplesToValueVector(datatable, KeyTypes, accumulated_keys);
	assert(valuevectors.size() == accumulated_keys.size());

	return valuevectors;
}

void TableDataOperator::insertValues(
		Table &_table) const
{
	// do nothing on empty keys
	if (accumulated_keys.empty())
		return;

	finalizeData();

	Table::Tuple_t &tuple = _table.getTuple();

	for (Table::TokenTypeMap_t::const_iterator iter = accumulatedData.begin();
			iter != accumulatedData.end();
			++iter) {
		tuple.replace( iter->first, iter->second);
	}
}


