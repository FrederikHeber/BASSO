/*
 * Database.hpp
 *
 *  Created on: Jun 12, 2014
 *      Author: heber
 */

#ifndef DATABASE_HPP_
#define DATABASE_HPP_

#include "BassoConfig.h"

#include <boost/shared_ptr.hpp>
#include <map>
#include <string>
#include <vector>

#include "Database/Table.hpp"

/** The Database class provides a single instance to contain all iteration-
 * related data. The data is eventually placed in a single-file (sqlite)
 * allowing for command-line parsing and data extraction for creating figures.
 */
class Database
{
public:
	//!> typedef for wrapping Database instance in shared_ptr
	typedef boost::shared_ptr<Database> Database_ptr_t;

	Database();
	~Database();

	void setDatabaseFile( const std::string &_filename)
	{ filename = _filename; DatabaseFileGiven = true; }

	/** Add a new table with name \a _name to this database.
	 *
	 * @param _name name of new table
	 * @return ref to the table, present one if it already existed
	 */
	Table& addTable(const std::string &_name);

	/** Removes a table with name \a _name from this database.
	 *
	 * @param _name name of new table
	 * @return true - table found and removed, false - else
	 */
	bool removeTable(const std::string &_name);

	/** Getter for a table by its name \a _name.
	 *
	 * @param _name name of table
	 * @return ref to the table
	 */
	Table& getTable(const std::string &_name);

	/** Getter for a table by its name \a _name.
	 *
	 * @param _name name of table
	 * @return const ref to the table
	 */
	const Table& getTableConst(const std::string &_name) const;

	/** Returns the number of tables present in the database.
	 *
	 * @return number of tables
	 */
	const size_t size() const
	{ return tables.size(); }

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

private:

	/** Internal function to write the database as an sqlite file.
	 *
	 * This takes the internal_database and writes it as a single
	 * table contained in a sqlite-compatible file.
	 *
	 * @return true - succesfully written, false - something went wrong
	 */
	bool writeSQLitefile() const;

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

	/** Writes a single table \a _table in the sqlite file.
	 *
	 * @param _table table to write to sqlite file
	 * @return true - everything ok, false - else
	 */
	bool writeTable(const Table &_table) const;

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

private:
	//!> states whether a database file name has been given or not.
	bool DatabaseFileGiven;

	//!> states whether we replace already present parameter tuples or just add new ones
	bool ReplacePresentParameterTuples;

	//!> contains the filename of the database for storing
	std::string filename;

	//!> static vector with all type names
	static std::vector<std::string> TypeNames;

	//!> typedef for the set of tables
	typedef std::map<std::string, Table> tables_t;

	//!> list of all tables of this database
	tables_t tables;

	//!> maximum number of keys per table
	const unsigned int MaxKeys;
};


#endif /* DATABASE_HPP_ */
