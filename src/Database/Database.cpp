/*
 * Database.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Database.hpp"

#include "boost/assign.hpp"

#include "Log/Logging.hpp"
#include "Poco/Data/Common.h"
#include "Poco/Data/SQLite/Connector.h"

using namespace Poco::Data;
using namespace boost::assign;

// static entities
std::vector<std::string> Database::TypeNames(Database::MAX_TYPES);


Database::Database() :
		DatabaseFileGiven(false)
{
	TypeNames[inttype] = "int";
	TypeNames[doubletype] = "double";
	TypeNames[valchartype] = "valchar(255)";
}

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

    if (!keys.empty()) {
		/// second, we need to gather a type for each key
		KeyType_t KeyTypes = getKeyToTypeMap(keys);
		assert(keys.size() == KeyTypes.size());

		/// third, we create the table with the given keys
		keys_t::const_iterator specialiter = keys.begin();
		const unsigned int MAXKEYS = 8;
		std::advance(specialiter, keys.size() > MAXKEYS ? MAXKEYS : keys.size());
		keys_t allowed_keys(keys.begin(), specialiter);
		Session ses("SQLite", filename.c_str());
		ses << "DROP TABLE IF EXISTS data", now;
		std::stringstream sql;
		{
			keys_t tempkeys(allowed_keys);
			sql << "CREATE TABLE data (";
			for (KeyType_t::const_iterator iter = KeyTypes.begin();
					iter != KeyTypes.end();) {
				if (tempkeys.find(iter->first) != tempkeys.end()) {
					sql << iter->first << " " << TypeNames[iter->second];
					tempkeys.erase(iter->first);
					if (!tempkeys.empty())
						sql << ",";
				}
				++iter;
			}
			sql << ")";
		}
		BOOST_LOG_TRIVIAL(trace)
			<< "SQL: " << sql.str();
		Statement stmt = ( ses << sql.str() );
		stmt.execute();

		/// then, convert the information in tuples in vector per type
		/// each having the same length
		std::vector<values_t> valuevector;
		{
			keys_t tempkeys(allowed_keys);
			for (KeyType_t::const_iterator keytypeiter = KeyTypes.begin();
					keytypeiter != KeyTypes.end(); ++keytypeiter) {
				if (tempkeys.find(keytypeiter->first) != tempkeys.end()) {
					switch(keytypeiter->second) {
					case inttype:
					{
						values_t values =
								getAllValuesPerType<int>(keytypeiter->first);
						valuevector.push_back(values);
						break;
					}
					case doubletype:
					{
						values_t values =
								getAllValuesPerType<double>(keytypeiter->first);
						valuevector.push_back(values);
						break;
					}
					case valchartype:
					{
						values_t values =
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
		}
		switch (valuevector.size()) {
		case 0:
			BOOST_LOG_TRIVIAL(warning)
				<< "Database contains no values";
			break;
		case 1:
			ses << "INSERT INTO data VALUES(?)", use(valuevector[0]), now;
			break;
		case 2:
			ses << "INSERT INTO data VALUES(?,?)", use(valuevector[0]), use(valuevector[1]), now;
			break;
		case 3:
			ses << "INSERT INTO data VALUES(?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), now;
			break;
		case 4:
			ses << "INSERT INTO data VALUES(?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), now;
			break;
		case 5:
			ses << "INSERT INTO data VALUES(?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), now;
			break;
		case 6:
			ses << "INSERT INTO data VALUES(?,?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), use(valuevector[5]), now;
			break;
		case 7:
			ses << "INSERT INTO data VALUES(?,?,?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), use(valuevector[5]), use(valuevector[6]), now;
			break;
		case 8:
			ses << "INSERT INTO data VALUES(?,?,?,?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), use(valuevector[5]), use(valuevector[6]), use(valuevector[7]), now;
			break;
		default:
			BOOST_LOG_TRIVIAL(error)
				<< "Cannot deal (yet) with tuple size" << valuevector.size();
			break;
		}

		/// finally, we write all information to the table
    } else {
    	BOOST_LOG_TRIVIAL(warning)
    			<< "The database is empty, not writing iteration-file.";
    }

    return true;
}

static
enum Database::types_t
getVariantsType(const boost::variant<int, double, std::string> &_variant)
{
	if (boost::get<const int>(&_variant) != NULL)
		return Database::inttype;
	else if (boost::get<const double>(&_variant) != NULL)
		return Database::doubletype;
	else if (boost::get<const std::string>(&_variant) != NULL)
		return Database::valchartype;
	else
		return Database::MAX_TYPES;
}

Database::KeyType_t Database::getKeyToTypeMap(
		const keys_t &_keys) const
{
	KeyType_t KeyType;

	for (internal_database_t::const_iterator iter = internal_database.begin();
			iter != internal_database.end(); ++iter) {
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
		for (Tuple_t::const_iterator keyiter = iter->begin();
				keyiter != iter->end(); ++keyiter) {
			// if key present in given ones
			if (_keys.count(keyiter->first)) {
				const enum types_t type = getVariantsType(keyiter->second);
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

Database::keys_t Database::getSetofUniqueKeys() const
{
	keys_t keys;

	for (internal_database_t::const_iterator iter = internal_database.begin();
			iter != internal_database.end(); ++iter) {
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

void Database::Tuple_t::replace(
		const std::string &_key,
		const typevariant_t &_value)
{
	std::pair<iterator, bool> inserter =
			insert( std::make_pair(_key, _value) );
	if (!inserter.second) {
		inserter.first->second = _value;
	}
}
