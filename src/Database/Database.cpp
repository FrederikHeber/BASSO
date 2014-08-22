/*
 * Database.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Database.hpp"

#include "boost/assign.hpp"
#include <sstream>

#include "Log/Logging.hpp"
#include "Poco/Data/Common.h"
#include "Poco/Data/SQLite/Connector.h"

#include "Table.hpp"

using namespace Poco::Data;
using namespace boost::assign;

// static entities
std::vector<std::string> Database::TypeNames(Database::MAX_TYPES);


Database::Database() :
		DatabaseFileGiven(false),
		MAXKEYS(12)
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

bool Database::writeSQLitefile()
{

	for (tables_t::const_iterator tableiter = tables.begin();
			tableiter != tables.end(); ++tableiter) {
		const Table &currenttable = tableiter->second;

		/// first, we need to find the set of unique keys
		const Table::keys_t keys = currenttable.getSetofUniqueKeys();

		if (!keys.empty()) {
			/// second, we need to gather a type for each key
			Table::KeyType_t KeyTypes = currenttable.getKeyToTypeMap(keys);
			assert(keys.size() == KeyTypes.size());

			/// third, we create the table with the given keys
			Table::keys_t::const_iterator specialiter = keys.begin();
			std::advance(specialiter, keys.size() > MAXKEYS ? MAXKEYS : keys.size());
			Table::keys_t allowed_keys(keys.begin(), specialiter);
			Session ses("SQLite", filename.c_str());
			// Don't drop table, we might want to accumulate multiple datasets
	//		ses << "DROP TABLE IF EXISTS " << currenttable.getName(), now;
			std::stringstream sql;
			{
				Table::keys_t tempkeys(allowed_keys);
				sql << "CREATE TABLE IF NOT EXISTS " << currenttable.getName() << " (";
				for (Table::KeyType_t::const_iterator iter = KeyTypes.begin();
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
			std::vector<Table::values_t> valuevector;
			{
				Table::keys_t tempkeys(allowed_keys);
				for (Table::KeyType_t::const_iterator keytypeiter = KeyTypes.begin();
						keytypeiter != KeyTypes.end(); ++keytypeiter) {
					if (tempkeys.find(keytypeiter->first) != tempkeys.end()) {
						switch(keytypeiter->second) {
						case inttype:
						{
							Table::values_t values =
									currenttable.getAllValuesPerType<int>(keytypeiter->first);
							valuevector.push_back(values);
							break;
						}
						case doubletype:
						{
							Table::values_t values =
									currenttable.getAllValuesPerType<double>(keytypeiter->first);
							valuevector.push_back(values);
							break;
						}
						case valchartype:
						{
							Table::values_t values =
									currenttable.getAllValuesPerType<std::string>(keytypeiter->first);
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
			ses << "BEGIN", now;
			switch (valuevector.size()) {
			case 0:
				BOOST_LOG_TRIVIAL(warning)
					<< "Database contains no values";
				break;
			case 1:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?)", use(valuevector[0]), now;
				break;
			case 2:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?,?)", use(valuevector[0]), use(valuevector[1]), now;
				break;
			case 3:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), now;
				break;
			case 4:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), now;
				break;
			case 5:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), now;
				break;
			case 6:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?,?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), use(valuevector[5]), now;
				break;
			case 7:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?,?,?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), use(valuevector[5]), use(valuevector[6]), now;
				break;
			case 8:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?,?,?,?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), use(valuevector[5]), use(valuevector[6]), use(valuevector[7]), now;
				break;
			case 9:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?,?,?,?,?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), use(valuevector[5]), use(valuevector[6]), use(valuevector[7]), use(valuevector[8]), now;
				break;
			case 10:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?,?,?,?,?,?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), use(valuevector[5]), use(valuevector[6]), use(valuevector[7]), use(valuevector[8]), use(valuevector[9]), now;
				break;
			case 11:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?,?,?,?,?,?,?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), use(valuevector[5]), use(valuevector[6]), use(valuevector[7]), use(valuevector[8]), use(valuevector[9]), use(valuevector[10]), now;
				break;
			case 12:
				ses << "INSERT INTO " << currenttable.getName() << " VALUES(?,?,?,?,?,?,?,?,?,?,?,?)", use(valuevector[0]), use(valuevector[1]), use(valuevector[2]), use(valuevector[3]), use(valuevector[4]), use(valuevector[5]), use(valuevector[6]), use(valuevector[7]), use(valuevector[8]), use(valuevector[9]), use(valuevector[10]), use(valuevector[11]), now;
				break;
			default:
				BOOST_LOG_TRIVIAL(error)
					<< "Cannot deal (yet) with tuple size" << valuevector.size();
				break;
			}
			ses << "END", now;

			/// finally, we write all information to the table
		} else {
			BOOST_LOG_TRIVIAL(warning)
					<< "The database is empty, not writing iteration-file.";
		}
	}

    return true;
}

Table& Database::addTable(const std::string &_name)
{
	// check whether table of such a name is not already present
	tables_t::iterator iter = tables.find(_name);
	if (iter == tables.end()) {
		tables.insert( std::make_pair(_name, Table(_name)) );
		iter = tables.find(_name);
	}
	return iter->second;
}

Table& Database::getTable(const std::string &_name)
{
	tables_t::iterator iter = tables.find(_name);
	assert( iter != tables.end() );
	return iter->second;
}

const Table& Database::getTableConst(const std::string &_name) const
{
	tables_t::const_iterator iter = tables.find(_name);
	assert( iter != tables.end() );
	return iter->second;
}

