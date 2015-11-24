/*
 * SQLDatabase.hpp
 *
 *  Created on: Nov 24, 2015
 *      Author: heber
 */

#ifndef SQLDATABASE_HPP_
#define SQLDATABASE_HPP_

#include "BassoConfig.h"

#include <string>
#include <vector>

#include "Database/TableDirectoryDatabase.hpp"

/** This class implements the Database interface as an sqlite-database.
 *
 */
class SQLDatabase : public TableDirectoryDatabase
{
public:
	SQLDatabase();
	virtual ~SQLDatabase();

	/** Clears the whole database.
	 *
	 */
	virtual void clear();

	void setDatabaseFile( const std::string &_filename)
	{ filename = _filename; DatabaseFileGiven = true; }

	/** Setter for whether tuples are replaced (with respect to their parameter
	 * part or we just add new ones regardless of same-parametered present ones.
	 *
	 * \warning Replacement can take quite some time for a larger database
	 * files (e.g. >10Mb ~ secs, >100Mb ~ mins).
	 *
	 * @param _flag true - replace, false - just add
	 */
	void setReplacePresentParameterTuples(const bool _flag)
	{ ReplacePresentParameterTuples = _flag; }

	/** Checks whether a tuple is present in the sqlite table without
	 * causing any warning if table is not uptodate.
	 *
	 * @param _table table to inspect
	 * @param _tuple tuple to look for
	 * @return true - tuple present, false - not present
	 */
	virtual bool isTuplePresentInTable(
			const Table &_table,
			const Table::Tuple_t &_tuple) const;

	/** Returns the rowid (sqlite specific primary key) for the requested
	 * \a _tuple in table \a _table.
	 *
	 * @param _table table to inspect
	 * @param _tuple tuple to look for
	 * @return rowid or (size_t)-1 of not present
	 */
	virtual size_t getIdOfTuplePresentInTable(
			const Table &_table,
			const Table::Tuple_t &_tuple) const;

	/** Writes a single table \a _table in the sqlite file.
	 *
	 * @param _table table to write to sqlite file
	 * @return true - everything ok, false - else
	 */
	virtual bool writeTable(const Table &_table) const;

	/** Execute a given string SQL command \a _command.
	 *
	 * \param _command SQL statement to execute
	 * \return true - everything ok, false - command failed
	 */
	virtual bool executeSQLStatement(const std::string &_command) const;

	/** Writes all tables in the database as an sqlite file.
	 *
	 * This takes the internal_database and writes it as a single
	 * table contained in a sqlite-compatible file.
	 *
	 * @return true - succesfully written, false - something went wrong
	 */
	virtual bool writeAllTables() const;

	/** Getter for whether a database file was given.
	 *
	 * This is important for whether tables can be written to file and thus
	 * whether a unique rowid id (for parameter sets) can be obtained or not.
	 *
	 * @return true - filename given, false - else
	 */
	bool isDatabaseFileGiven() const
	{ return DatabaseFileGiven; }

private:

	/** Adds the table to the sqlite file if it not already exists.
	 *
	 * @param _table table to create
	 * @param _KeyTypes map with type for each key (column name)
	 * @param _allowed_keys vector of all allowed keys (filter for column names)
	 * @return true - everything ok, false - else
	 */
	bool createTableIfNotExists(
			const Table &_table,
			const Table::KeyType_t &_KeyTypes,
			const Table::keys_t &_allowed_keys) const;

	/** Deletes all present tuples in given \a _table the sqlite file.
	 *
	 * @param _table table whose present tuples to remove
	 * @param _KeyTypes map with type for each key (column name)
	 * @return true - everything ok, false - else
	 */
	bool deletePresentTuplesinTable(
			const Table &_table,
			const Table::KeyType_t &_KeyTypes) const;

	/** Reads a single table \a _table from the sqlite file.
	 *
	 * \note We can only parse columns that are already present in \a
	 * _table. This is basically a limitation of the poco library we
	 * use to get sqlite functionality. Sqlite knows about
	 * \code
	 * 	PRAGMA table_info(table_name);
	 * \endcode
	 * but we cannot squeeze the information into an STL container; it
	 * will just go into void. Hence, the tables is deleted (if non-
	 * empty) and filled with all tuples from columns that were present
	 * already before.
	 *
	 * @param _table table to read
	 * @return true - everything ok, false - else
	 */
	bool readTable(Table &_table);

	/** Updates the values in the sqlite table with the same named
	 * table \a _table.
	 *
	 * @param _table table to take values from for update
	 * @param _KeyTypes column name to type map
	 * @param _allowed_keys set of keys (column names) to modify
	 * @return true - everything ok, false - else
	 */
	bool updateTable(
			const Table &_table,
			const Table::KeyType_t &_KeyTypes,
			const Table::keys_t &_allowed_keys) const;

	std::string printTupleWhereStatement(
			const Table::Tuple_t &_tuple,
			const Table::KeyType_t &_KeyTypes) const;

	bool hasTupleOneParameter(
			const Table::Tuple_t &_tuple) const;

private:

	//!> states whether a database file name has been given or not.
	bool DatabaseFileGiven;

	//!> states whether we replace already present parameter tuples or just add new ones
	bool ReplacePresentParameterTuples;

	//!> contains the filename of the database for storing
	std::string filename;

	//!> static vector with all type names
	static std::vector<std::string> TypeNames;

	//!> maximum number of keys per table
	const unsigned int MaxKeys;
};


#endif /* SQLDATABASE_HPP_ */
