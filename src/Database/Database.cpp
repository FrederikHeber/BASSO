/*
 * Database.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Database.hpp"

#include "Poco/Data/Common.h"
#include "Poco/Data/SQLite/Connector.h"

using namespace Poco::Data;

Database::Database() :
		DatabaseFileGiven(false)
{}

Database::~Database()
{
	// open connection
	SQLite::Connector::registerConnector();

	// write information
	if (DatabaseFileGiven) {
		writeSQLitefile();
	}

    // quit
	SQLite::Connector::unregisterConnector();
}

void Database::addTuple(const Tuple_t &_tuple)
{
	internal_database.insert(_tuple);
}

bool Database::writeSQLitefile()
{
	/// first, we need to find the set of unique keys
	const keys_t keys = getSetofUniqueKeys();

	/// second, we need to gather a type for each key
	KeyType_t KeyTypes = getKeyToTypeMap(keys);
    assert(keys.size() == KeyTypes.size());

	/// third, we create the table with the given keys
    Session ses("SQLite", filename.c_str());
    ses << "CREATE TABLE data (";
    for (KeyType_t::const_iterator iter = KeyTypes.begin();
    		iter != KeyTypes.end();) {
		ses << iter->first << " " << iter->second;
		++iter;
		if (iter != KeyTypes.end())
			ses << ",";
    }
    ses << ");";

    /// then, convert the information in tuples in vector per type
    /// each having the same length
    values_t values =
    		getAllValuesPerType(keys_t(keys.begin(), ++keys.begin()));

	/// finally, we write all information to the table
    ses << "INSERT INTO data VALUES(:" << *keys.begin()
    		<< ")", use(values), now;
}

Database::values_t
Database::getAllValuesPerType(const keys_t &_keys) const
{
	values_t values;

	for (keys_t::const_iterator keyiter = _keys.begin();
			keyiter != _keys.end(); ++keyiter) {
		for (internal_database_t::const_iterator iter = internal_database.begin();
				iter != internal_database.end(); ++iter) {
			const Tuple_t &currrent = *iter;
			Tuple_t::const_iterator finditer = iter->find(*keyiter);
			// if key present in given ones
			if (finditer != iter->end()) {
				values.push_back(boost::get<double>(finditer->second));
			} else
				values.push_back(0.);
		}
	}

	return values;
}

const std::string&
getVariantsType(const boost::variant<int, double, std::string> &_variant)
{
	static const std::string intname = "int";
	static const std::string doublename = "double";
	static const std::string stringname = "varchar(255)";
	static const std::string unknown = "unknown";
	if (boost::get<const int>(&_variant) != NULL)
		return intname;
	else if (boost::get<const double>(&_variant) != NULL)
		return doublename;
	else if (boost::get<const std::string>(&_variant) != NULL)
		return stringname;
	else
		return unknown;
}

Database::KeyType_t Database::getKeyToTypeMap(
		const keys_t &_keys) const
{
	typedef std::map< std::string, std::string > KeyType_t;
	KeyType_t KeyType;

	for (internal_database_t::const_iterator iter = internal_database.begin();
			iter != internal_database.end(); ++iter) {
		const Tuple_t &currrent = *iter;
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

bool Database::checkDatabaseSanity(const keys_t &_keys) const
{
	KeyType_t KeyType;
	bool status = true;

	for (internal_database_t::const_iterator iter = internal_database.begin();
			iter != internal_database.end(); ++iter) {
		const Tuple_t &currrent = *iter;
		for (Tuple_t::const_iterator keyiter = iter->begin();
				keyiter != iter->end(); ++keyiter) {
			// if key present in given ones
			if (_keys.count(keyiter->first)) {
				const std::string &NameOfType = getVariantsType(keyiter->second);
				std::pair< KeyType_t::iterator, bool > inserter =
						KeyType.insert( std::make_pair(
								keyiter->first,
								NameOfType)
						);
				// if value present check against map
				if (inserter.second == false)
					status &= inserter.first->second == NameOfType;
			}
		}
	}
	return status;
}

Database::keys_t Database::getSetofUniqueKeys() const
{
	keys_t keys;

	for (internal_database_t::const_iterator iter = internal_database.begin();
			iter != internal_database.end(); ++iter) {
		const Tuple_t &currrent = *iter;
		for (Tuple_t::const_iterator keyiter = iter->begin();
				keyiter != iter->end(); ++keyiter) {
			keys.insert(keyiter->first);
		}
	}

	return keys;
}

bool Database::Tuple_t::operator<(const Tuple_t &_a) const
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
