/*
 * Database.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Database.hpp"

#include <boost/assign.hpp>
#include <sstream>

#include "Log/Logging.hpp"
#include "Poco/Data/Common.h"
#include "Poco/Data/SQLite/Connector.h"

#include "Table.hpp"

using namespace Poco::Data;
using namespace boost::assign;

// static entities
std::vector<std::string> Database::TypeNames(Database::MAX_TYPES);


#define BASSO_MAXKEYS 16

Database::Database() :
		DatabaseFileGiven(false),
		ReplacePresentParameterTuples(false),
		MaxKeys(BASSO_MAXKEYS)
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
			std::advance(specialiter, keys.size() > MaxKeys ? MaxKeys : keys.size());
			if (keys.size() > MaxKeys)
				BOOST_LOG_TRIVIAL(error)
						<< "Truncating database output by " << (keys.size() - MaxKeys)
						<< " columns as too many columns are to be stored.";
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

			if (ReplacePresentParameterTuples) {
				/// convert all tuples associated with parameters to pairs
				/// containing column name and associated value, this is
				/// required for specifying which present rows to remove
				ses << "BEGIN", now;
				for (Table::internal_table_t::const_iterator tupleiter = tableiter->second.internal_table.begin();
						tupleiter != tableiter->second.internal_table.end(); ++tupleiter) {
					std::stringstream sql_delete;
					sql_delete << "DELETE FROM " << currenttable.getName()
						<< " WHERE ";
					bool AtLeastOneParameter = false;
					for (Table::Tuple_t::const_iterator paramiter = tupleiter->begin();
							paramiter != tupleiter->end(); ++paramiter) {
						if (tupleiter->isParameter(paramiter->first)) {
							if (AtLeastOneParameter)
								sql_delete << " AND ";
							sql_delete << paramiter->first << "=";
							Table::KeyType_t::const_iterator keytypeiter =
									KeyTypes.find(paramiter->first);
							assert( keytypeiter != KeyTypes.end() );
							switch(keytypeiter->second) {
							case valchartype:
								sql_delete << "'" << paramiter->second << "'";
								break;
							default:
								sql_delete << paramiter->second;
								break;
							}
							AtLeastOneParameter = true;
						}
					}
					sql_delete << ";";
					// only execute if WHERE statement makes sense
					if (AtLeastOneParameter) {
						ses << sql_delete.str(), now;
						BOOST_LOG_TRIVIAL(trace)
							<< "SQL: " << sql_delete.str();
					}
				}
				ses << "END", now;
			}

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

			// add all new ones
			// note that due to use(..) this cannot be debugged into
			// a printable string
			{
				switch (valuevector.size()) {
				case 0:
					BOOST_LOG_TRIVIAL(warning)
						<< "Database contains no values";
					break;
				// use preprocessor magic to create the range of cases
#include <boost/preprocessor/iteration/local.hpp>
#include "Database_impl.hpp"
#define BASSO_ARGUMENTLIST (ses)(currenttable.getName())(valuevector)
#define BOOST_PP_LOCAL_MACRO(n) CasePrinter(~, n, BASSO_ARGUMENTLIST)
#define BOOST_PP_LOCAL_LIMITS (1, BASSO_MAXKEYS)
#include BOOST_PP_LOCAL_ITERATE()
#undef BASSO_ARGUMENTLIST
#include "Database_undef.hpp"
				default:
					BOOST_LOG_TRIVIAL(error)
						<< "Cannot deal (yet) with tuple size" << valuevector.size();
					break;
				}
			}

			ses << "END", now;

			/// finally, we write all information to the table
		} else {
			BOOST_LOG_TRIVIAL(warning)
					<< "The table " << tableiter->first
					<< " is empty, not writing to iteration-file.";
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

