/*
 * Table.cpp
 *
 *  Created on: Aug 22, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Database.hpp"
#include "Table.hpp"

#include "Log/Logging.hpp"


void Table::addTuple(const Tuple_t &_tuple)
{
	uptodate = false;
	internal_table.insert(_tuple);
}

static
enum Database_types::types_t
getVariantsType(const boost::variant<int, double, std::string> &_variant)
{
	if (boost::get<const int>(&_variant) != NULL)
		return Database_types::inttype;
	else if (boost::get<const double>(&_variant) != NULL)
		return Database_types::doubletype;
	else if (boost::get<const std::string>(&_variant) != NULL)
		return Database_types::valchartype;
	else
		return Database_types::MAX_TYPES;
}

Table::KeyType_t Table::getKeyToTypeMap(
		const keys_t &_keys) const
{
	KeyType_t KeyType;

	for (internal_table_t::const_iterator iter = internal_table.begin();
			iter != internal_table.end(); ++iter) {
		for (Tuple_t::const_iterator keyiter = iter->begin();
				keyiter != iter->end(); ++keyiter) {
			// if key present in given ones
			if (_keys.count(keyiter->first)) {
				// we ignore if the value is already present
				KeyType.insert(
						std::make_pair(
								keyiter->first,
								getVariantsType(keyiter->second)
								)
					);
			}
		}
	}

	return KeyType;
}

bool Table::checkTableSanity(const keys_t &_keys) const
{
	KeyType_t KeyType;
	bool status = true;

	for (internal_table_t::const_iterator iter = internal_table.begin();
			iter != internal_table.end(); ++iter) {
		for (Tuple_t::const_iterator keyiter = iter->begin();
				keyiter != iter->end(); ++keyiter) {
			// if key present in given ones
			if (_keys.count(keyiter->first)) {
				const enum Database_types::types_t type =
						getVariantsType(keyiter->second);
				std::pair< KeyType_t::iterator, bool > inserter =
						KeyType.insert( std::make_pair(
								keyiter->first,
								type)
						);
				// if value present check against map
				if (inserter.second == false)
					status &= inserter.first->second == type;
			}
		}
	}
	return status;
}

Table::keys_t Table::getSetofUniqueKeys() const
{
	keys_t keys;

	for (internal_table_t::const_iterator iter = internal_table.begin();
			iter != internal_table.end(); ++iter) {
		for (Tuple_t::const_iterator keyiter = iter->begin();
				keyiter != iter->end(); ++keyiter) {
			keys.insert(keyiter->first);
		}
	}

	return keys;
}

std::vector<Table::values_t> Table::convertTuplesToValueVector(
		const KeyType_t &_KeyTypes,
		const keys_t &_allowed_keys) const
{
	std::vector<Table::values_t> valuevector;
	valuevector.reserve(_allowed_keys.size());

	Table::keys_t tempkeys(_allowed_keys);
	for (Table::KeyType_t::const_iterator keytypeiter = _KeyTypes.begin();
			keytypeiter != _KeyTypes.end(); ++keytypeiter) {
		if (tempkeys.find(keytypeiter->first) != tempkeys.end()) {
			switch(keytypeiter->second) {
			case Database_types::inttype:
			{
				Table::values_t values =
						getAllValuesPerType<int>(keytypeiter->first);
				valuevector.push_back(values);
				break;
			}
			case Database_types::doubletype:
			{
				Table::values_t values =
						getAllValuesPerType<double>(keytypeiter->first);
				valuevector.push_back(values);
				break;
			}
			case Database_types::valchartype:
			{
				Table::values_t values =
						getAllValuesPerType<std::string>(keytypeiter->first);
				valuevector.push_back(values);
				break;
			}
			default:
				BOOST_LOG_TRIVIAL(error)
					<< "Unknown type for key " << keytypeiter->first;
				break;
			}
			tempkeys.erase(keytypeiter->first);
		}
	}

	return valuevector;
}

void Table::clear()
{
	reset();
}

void Table::reset()
{
	internal_table.clear();
	internal_tuple.clear();
	uptodate = true;
}
