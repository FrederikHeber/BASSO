/*
 * Database.cpp
 *
 *  Created on: Jun 12, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Database.hpp"

#include <boost/assign.hpp>
//#include <iterator>
#include <sstream>

#include "Log/Logging.hpp"
#include "Poco/Data/Common.h"
#include "Poco/Data/SQLite/Connector.h"

#include "Table.hpp"

using namespace Poco::Data;
using namespace boost::assign;

// static entities
std::vector<std::string> Database::TypeNames(Database_types::MAX_TYPES);

#define BASSO_MAXKEYS 18

Database::Database() :
		DatabaseFileGiven(false),
		ReplacePresentParameterTuples(false),
		MaxKeys(BASSO_MAXKEYS)
{
	TypeNames[Database_types::inttype] = "int";
	TypeNames[Database_types::doubletype] = "double";
	TypeNames[Database_types::valchartype] = "valchar(255)";
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

bool Database::writeSQLitefile() const
{

	for (tables_t::const_iterator tableiter = tables.begin();
			tableiter != tables.end(); ++tableiter) {
		const Table &currenttable = *tableiter->second;
		writeTable(currenttable);
	}

    return true;
}

bool Database::createTableIfNotExists(
		const Table &_table,
		const Table::KeyType_t &_KeyTypes,
		const Table::keys_t &_allowed_keys) const
{
	Session ses("SQLite", filename.c_str());
	// Don't drop table, we might want to accumulate multiple datasets
//		ses << "DROP TABLE IF EXISTS " << _table.getName(), now;
	std::stringstream sql;
	{
		Table::keys_t tempkeys(_allowed_keys);
		sql << "CREATE TABLE IF NOT EXISTS " << _table.getName() << " (";
		for (Table::KeyType_t::const_iterator iter = _KeyTypes.begin();
				iter != _KeyTypes.end();) {
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

	return true;
}

bool Database::deletePresentTuplesinTable(
		const Table &_table,
		const Table::KeyType_t &_KeyTypes) const
{
	/// convert all tuples associated with parameters to pairs
	/// containing column name and associated value, this is
	/// required for specifying which present rows to remove
	Session ses("SQLite", filename.c_str());
	ses << "BEGIN", now;
	for (Table::internal_table_t::const_iterator tupleiter = _table.internal_table.begin();
			tupleiter != _table.internal_table.end(); ++tupleiter) {
		std::stringstream sql_delete;
		sql_delete << "DELETE FROM " << _table.getName()
			<< " WHERE ";
		bool AtLeastOneParameter = false;
		for (Table::Tuple_t::const_iterator paramiter = tupleiter->begin();
				paramiter != tupleiter->end(); ++paramiter) {
			if (tupleiter->isParameter(paramiter->first)) {
				if (AtLeastOneParameter)
					sql_delete << " AND ";
				sql_delete << paramiter->first << "=";
				Table::KeyType_t::const_iterator keytypeiter =
						_KeyTypes.find(paramiter->first);
				assert( keytypeiter != _KeyTypes.end() );
				switch(keytypeiter->second) {
				case Database_types::valchartype:
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

	return true;
}

bool Database::updateTable(
		const Table &_table,
		const Table::KeyType_t &_KeyTypes,
		const Table::keys_t &_allowed_keys) const
{
	// first, convert the information in tuples in vector per type
	// each having the same length
	std::vector<Table::values_t> valuevector =
			_table.convertTuplesToValueVector(_KeyTypes, _allowed_keys);

	Session ses("SQLite", filename.c_str());
	ses << "BEGIN", now;

	// then, add all new ones
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
#define BASSO_ARGUMENTLIST (ses)(_table.getName())(valuevector)
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

	// finally, we write all information to the table, i.e. end transaction
	ses << "END", now;

	// set uptodate flag
	_table.uptodate = true;

	return true;
}

bool Database::writeTable(const Table &_table) const
{
	// skip if table is uptodate
	if (_table.uptodate)
		return true;

	/// first, we need to find the set of unique keys
	const Table::keys_t keys = _table.getSetofUniqueKeys();

	if (!keys.empty()) {
		/// second, we need to gather a type for each key
		Table::KeyType_t KeyTypes = _table.getKeyToTypeMap(keys);
		assert(keys.size() == KeyTypes.size());

		/// third, we create the table with the given keys
		Table::keys_t::const_iterator specialiter = keys.begin();
		std::advance(specialiter, keys.size() > MaxKeys ? MaxKeys : keys.size());
		if (keys.size() > MaxKeys)
			BOOST_LOG_TRIVIAL(error)
					<< "Truncating database output by " << (keys.size() - MaxKeys)
					<< " columns as too many columns are to be stored.";
		Table::keys_t allowed_keys(keys.begin(), specialiter);
		createTableIfNotExists(_table, KeyTypes, allowed_keys);

		if (ReplacePresentParameterTuples)
			deletePresentTuplesinTable(_table, KeyTypes);

		/// Combine the following update into single transaction
		updateTable(_table, KeyTypes, allowed_keys);

	} else {
		BOOST_LOG_TRIVIAL(warning)
				<< "The table " << _table.getName()
				<< " is empty, not writing to iteration-file.";
	}

	return true;
}

bool Database::readTable(
		Table& _table)
{
	/// check if table is uptodate
	if (!_table.isUptodate()) {
		BOOST_LOG_TRIVIAL(error)
				<< "Table " << _table.getName() << " has not been written.";
		return false;
	}

	// first, we need to keep the set of unique keys
	const Table::keys_t keys = _table.getSetofUniqueKeys();
	if (!keys.empty()) {
		// second, we need to keep the type for each key fixed
		Table::KeyType_t KeyTypes = _table.getKeyToTypeMap(keys);
		assert(keys.size() == KeyTypes.size());

		// clear table if not empty
		if (!_table.empty())
			_table.clear();

		// parse tuples from sqlite database
		Session ses("SQLite", filename.c_str());
		typedef std::map<std::string, Table::values_t > KeyValueMap_t;
		KeyValueMap_t KeyValueMap;
		for (Table::keys_t::const_iterator iter = keys.begin();
				iter != keys.end(); ++iter) {
			Table::values_t valuevector;
			ses << "SELECT " << *iter
				<< " FROM " << _table.getName(), into(valuevector), now;
//			BOOST_LOG_TRIVIAL(trace)
//				<< "SQL: " << "SELECT " << *iter
//				<< " FROM " << _table.getName();
//			std::stringstream valuestream;
//			std::copy(valuevector.begin(),valuevector.end(),
//					std::ostream_iterator<std::string>(valuestream, ","));
//			BOOST_LOG_TRIVIAL(trace)
//				<< "RESULT: #" << valuevector.size() << " elements ("
//				<< valuestream.str() << ")";
			KeyValueMap.insert( std::make_pair(*iter, valuevector) );
		}

		if (KeyValueMap.empty()) {
			BOOST_LOG_TRIVIAL(error)
					<< "Table " << _table.getName() << " is empty in sqlite file.";
			return false;
		}

		// check that each entry's vector in map has same size
		KeyValueMap_t::const_iterator checkiter = KeyValueMap.begin();
		const size_t firstsize = checkiter->second.size();
		assert(firstsize != 0);
		for (++checkiter; checkiter != KeyValueMap.end(); ++checkiter)
			if (firstsize != checkiter->second.size()) {
				BOOST_LOG_TRIVIAL(error)
						<< "Table " << _table.getName() << "'s column "
						<< checkiter->first << " has different number of rows ("
						<< checkiter->second.size() << ") than first column "
						<< KeyValueMap.begin()->first << "( "
						<< firstsize << ")";
				return false;
			}

		// initialize the tuple from first row
		Table::Tuple_t tuple;
		for (KeyValueMap_t::const_iterator iter = KeyValueMap.begin();
				iter != KeyValueMap.end(); ++iter) {
			assert( KeyTypes.count(iter->first) );
			std::stringstream valuestream(iter->second[firstsize-1]);
			switch (KeyTypes[iter->first]) {
			case Database_types::inttype:
			{
				int tempval;
				valuestream >> tempval;
				tuple.insert(
						std::make_pair(
								iter->first,
								tempval),
						Table::Data);
				break;
			}
			case Database_types::doubletype:
			{
				double tempval;
				valuestream >> tempval;
				tuple.insert(
						std::make_pair(
								iter->first,
								tempval),
						Table::Data);
				break;
			}
			case Database_types::valchartype:
			{
				tuple.insert(
						std::make_pair(
								iter->first,
								iter->second[firstsize-1]),
						Table::Data);
				break;
			}
			default:
				BOOST_LOG_TRIVIAL(error)
					<< "Key " << iter->first << " has unknown type.";
				return false;
			}
		}
		_table.addTuple(tuple);
		// and insert each tuple
		for (int row = firstsize-2; row >= 0; --row) {
			for (KeyValueMap_t::const_iterator iter = KeyValueMap.begin();
					iter != KeyValueMap.end(); ++iter) {
				assert( KeyTypes.count(iter->first) );
				std::stringstream valuestream(iter->second[row]);
				switch (KeyTypes[iter->first]) {
				case Database_types::inttype:
				{
					int tempval;
					valuestream >> tempval;
					tuple.replace(iter->first, tempval);
					break;
				}
				case Database_types::doubletype:
				{
					double tempval;
					valuestream >> tempval;
					tuple.replace(iter->first, tempval);
					break;
				}
				case Database_types::valchartype:
				{
					tuple.replace(iter->first, iter->second[row]);
					break;
				}
				default:
					BOOST_LOG_TRIVIAL(error)
						<< "Key " << iter->first << " has unknown type.";
					return false;
				}

			}
			_table.addTuple(tuple);
		}
	} else {
		BOOST_LOG_TRIVIAL(error)
				<< "Table " << _table.getName() << " has no keys in sqlite file.";
		return false;
	}

	return true;
}

Table& Database::addTable(const std::string &_name)
{
	// check whether table of such a name is not already present
	tables_t::iterator iter = tables.find(_name);
	if (iter == tables.end()) {
		Table::ptr insert_table(new Table(_name));
		tables.insert( std::make_pair(_name, insert_table) );
		iter = tables.find(_name);
	}
	return *iter->second;
}

bool Database::removeTable(const std::string &_name)
{
	tables_t::iterator iter = tables.find(_name);
	const bool status = iter != tables.end();
	if (status)
		tables.erase(iter);
	return status;
}

Table& Database::getTable(const std::string &_name)
{
	tables_t::iterator iter = tables.find(_name);
	assert( iter != tables.end() );
	return *iter->second;
}

const Table& Database::getTableConst(const std::string &_name) const
{
	tables_t::const_iterator iter = tables.find(_name);
	assert( iter != tables.end() );
	return *iter->second;
}

