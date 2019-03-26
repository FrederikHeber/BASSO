/*
 * TableDataAccumulator.cpp
 *
 *  Created on: Nov 24, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "TableDataAccumulator.hpp"

#include <algorithm>

#include "Database/Table.hpp"
#include "Log/Logging.hpp"
#include "MatrixFactorizer/Database/MaxMinAverage.hpp"

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

template <class T>
void insertValue(
		Table::TokenTypeMap_t &_accumulatedData,
		const std::string &_keyname,
		const Table::any_values_t &_valuevector
		)
{
	MaxMinAverage<T> Bregman(
			convertAnyVectorTo<T>(_valuevector));
	_accumulatedData["min_"+_keyname] = Bregman.minimum;
	_accumulatedData["max_"+_keyname] = Bregman.maximum;
	_accumulatedData["avg_"+_keyname] = Bregman.average;
	_accumulatedData["var_"+_keyname] = Bregman.variance;
}

void TableDataAccumulator::extractData() const
{
	// do nothing on empty keys
	if (accumulated_keys.empty())
		return;

	const valuevectors_t valuevectors = extractValues();
	for (Table::keys_t::const_iterator keyiter = accumulated_keys.begin();
			keyiter != accumulated_keys.end();
			++keyiter) {
		const std::string &keyname = *keyiter;
		const size_t iterindex =
				std::distance(accumulated_keys.begin(), keyiter);
		LOG(trace, "There are " << valuevectors[iterindex].second.size() << " values for " << keyname);
		assert(valuevectors[iterindex].second.size()==1);
		const AccumulatedValues::iterator iter =
				accumulatedValues.find(keyname);
		if (iter != accumulatedValues.end())
			iter->second.second.push_back(valuevectors[iterindex].second[0]);
		else
			accumulatedValues[keyname] =
					std::make_pair(
							valuevectors[iterindex].first,
							Table::any_values_t(1, valuevectors[iterindex].second[0])
					);
	}
}

void TableDataAccumulator::finalizeData() const
{
	for (AccumulatedValues::const_iterator valueiter = accumulatedValues.begin();
			valueiter != accumulatedValues.end(); ++valueiter) {
		const std::string &keyname = valueiter->first;
		const Database_types::types_t &type = valueiter->second.first;
		const Table::any_values_t &values = valueiter->second.second;
		switch(type) {
			case Database_types::inttype:
			{
				insertValue<int>(accumulatedData, keyname, values);
				break;
			}
			case Database_types::doubletype:
			{
				insertValue<double>(accumulatedData, keyname, values);
				break;
			}
			case Database_types::valchartype:
			default:
				LOG(error, "Unknown type for key " << keyname);
				break;
		}
	}
}

template <class T>
void insertDefaultValue(
		Table::Tuple_t &_tuple,
		const std::string &_keyname
		)
{
	DefaultValue<T> defaultvalue;
	_tuple.insert( std::make_pair("min_"+_keyname, defaultvalue.get()), Table::Data);
	_tuple.insert( std::make_pair("max_"+_keyname, defaultvalue.get()), Table::Data);
	_tuple.insert( std::make_pair("avg_"+_keyname, defaultvalue.get()), Table::Data);
	_tuple.insert( std::make_pair("var_"+_keyname, defaultvalue.get()), Table::Data);
}

void TableDataAccumulator::prepareTableForAccumulatedValues(
		Table &_table,
		const Database_types::types_t &_type,
		const std::string &_keyname)
{
	Table::Tuple_t &tuple = _table.getTuple();

	switch(_type) {
		case Database_types::inttype:
		{
			insertDefaultValue<int>(tuple, _keyname);
			break;
		}
		case Database_types::doubletype:
		{
			insertDefaultValue<double>(tuple, _keyname);
			break;
		}
		case Database_types::valchartype:
		default:
			LOG(error, "Unknown or unhandleable type for key " << _keyname);
			break;
	}
}

void TableDataAccumulator::addAccumulatedValues(
		const AccumulatedValues &_values)
{
	if (_values.size() != 0)
		accumulatedValues.insert(_values.begin(),_values.end());
}

