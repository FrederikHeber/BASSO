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
 *
 * \section Outline of the class and usage outline
 *
 * The Database is used by adding tables via addTable(), which immediately
 * returns a reference to the instantiated table. The reference of the new
 * table can retrieved any time via getTable() via its unique table name.
 * After the table has been filled with tuples, writeSQLitefile() writes
 * the whole set of known tables to disc. For updating just a single
 * table writeTable() can be used.
 *
 * \section A brief walk in SQL land sql-walkthru
 *
 * In general, SQL databases have multiple tables. They are linked via keys
 * where each row in a table is identified by a unique so-called \b primary
 * key. These keys may also be stored in other tables as so-called \b foreign
 * keys. Via matching primary and foreign keys two (or more) tables can
 * be linked in a so-called \b join to create a larger, combined table.
 * Last but not least, a so-called \b view can be thought of like a filter
 * that for example performs the stored join magically in the background,
 * while the user in the foreground simply operates with a normal table. In
 * this case the larger, combined table directly.
 *
 * \section Central idea idea
 *
 * The central idea for the database in a laboratory code is that you have
 * (at least) two tables: one table stores the parameter set and another
 * all the data -- e.g. iteration counts, residual information, ... for an
 * iterative method. The parameter set is identified via its (row)id and
 * this id is stored along with all associated (iterative) data in the
 * second table. Multiple runs with varying parameter sets differ by this
 * foreign key in the data table
 *
 */
class Database
{
public:
	//!> typedef for wrapping Database instance in shared_ptr
	typedef boost::shared_ptr<Database> Database_ptr_t;

	Database();
	virtual ~Database();

	void setDatabaseFile( const std::string &_filename)
	{ filename = _filename; DatabaseFileGiven = true; }

	/** Add a new table with name \a _name to this database.
	 *
	 * @param _name name of new table
	 * @return ref to the table, present one if it already existed
	 */
	virtual Table& addTable(const std::string &_name);

	/** Removes a table with name \a _name from this database.
	 *
	 * @param _name name of new table
	 * @return true - table found and removed, false - else
	 */
	virtual bool removeTable(const std::string &_name);

	/** Removes contents of a table with name \a _name from this database.
	 *
	 * @param _name name of new table
	 * @return true - table found and contents removed, false - else
	 */
	bool clearTable(const std::string &_name);

	/** Getter for a table by its name \a _name.
	 *
	 * @param _name name of table
	 * @return ref to the table
	 */
	virtual Table& getTable(const std::string &_name);

	/** Getter for a table by its name \a _name.
	 *
	 * @param _name name of table
	 * @return const ref to the table
	 */
	virtual const Table& getTableConst(const std::string &_name) const;

	/** Returns the number of tables present in the database.
	 *
	 * @return number of tables
	 */
	virtual const size_t size() const
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

	//!> typedef for the set of tables
	typedef std::map<std::string, Table::ptr> tables_t;

	//!> list of all tables of this database
	tables_t tables;

	//!> maximum number of keys per table
	const unsigned int MaxKeys;
};


#endif /* DATABASE_HPP_ */
