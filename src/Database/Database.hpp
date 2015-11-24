/*
 * Database.hpp
 *
 *  Created on: Jun 12, 2014
 *      Author: heber
 */

#ifndef DATABASE_HPP_
#define DATABASE_HPP_

#include "BassoConfig.h"

#include <string>

#include <boost/shared_ptr.hpp>

#include "Database/Table.hpp"

/** The Database class provides the interface for a single instance to contain
 * all iteration-related data. The data is eventually placed in a single-file
 * (sqlite) allowing for command-line parsing and data extraction for creating
 * figures.
 *
 * \section outline Outline of the class and usage outline
 *
 * The Database is used by adding tables via addTable(), which immediately
 * returns a reference to the instantiated table. The reference of the new
 * table can retrieved any time via getTable() via its unique table name.
 * After the table has been filled with tuples, writeSQLitefile() writes
 * the whole set of known tables to disc. For updating just a single
 * table writeTable() can be used.
 *
 * \section sql-walk A brief walk in SQL land sql-walkthru
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
 * \section central-idea Central idea
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

	/** Cstor of class Database.
	 *
	 */
	Database()
	{}

	/** Dstor of class Database.
	 *
	 */
	virtual ~Database()
	{}

	/** Set the filename where the database is stored.
	 *
	 * @param _filename file name to write database to
	 */
	virtual void setDatabaseFile( const std::string &_filename) = 0;

	/** Add a new table with name \a _name to this database.
	 *
	 * @param _name name of new table
	 * @return ref to the table, present one if it already existed
	 */
	virtual Table& addTable(const std::string &_name) = 0;

	/** Removes a table with name \a _name from this database.
	 *
	 * @param _name name of new table
	 * @return true - table found and removed, false - else
	 */
	virtual bool removeTable(const std::string &_name) = 0;

	/** Removes contents of a table with name \a _name from this database.
	 *
	 * @param _name name of new table
	 * @return true - table found and contents removed, false - else
	 */
	virtual bool clearTable(const std::string &_name) = 0;

	/** Getter for a table by its name \a _name.
	 *
	 * @param _name name of table
	 * @return ref to the table
	 */
	virtual Table& getTable(const std::string &_name) = 0;

	/** Getter for a table by its name \a _name.
	 *
	 * @param _name name of table
	 * @return const ref to the table
	 */
	virtual const Table& getTableConst(const std::string &_name) const = 0;

	/** Returns the number of tables present in the database.
	 *
	 * @return number of tables
	 */
	virtual const size_t size() const = 0;

	/** Clears the whole database.
	 *
	 */
	virtual void clear() = 0;

	/** Setter for whether tuples are replaced (with respect to their parameter
	 * part or we just add new ones regardless of same-parametered present ones.
	 *
	 * \warning Replacement can take quite some time for a larger database
	 * files (e.g. >10Mb ~ secs, >100Mb ~ mins).
	 *
	 * @param _flag true - replace, false - just add
	 */
	virtual void setReplacePresentParameterTuples(const bool _flag) = 0;

	/** Checks whether a tuple is present in the sqlite table without
	 * causing any warning if table is not uptodate.
	 *
	 * @param _table table to inspect
	 * @param _tuple tuple to look for
	 * @return true - tuple present, false - not present
	 */
	virtual bool isTuplePresentInTable(
			const Table &_table,
			const Table::Tuple_t &_tuple) const = 0;

	/** Returns the rowid (sqlite specific primary key) for the requested
	 * \a _tuple in table \a _table.
	 *
	 * @param _table table to inspect
	 * @param _tuple tuple to look for
	 * @return rowid or (size_t)-1 of not present
	 */
	virtual size_t getIdOfTuplePresentInTable(
			const Table &_table,
			const Table::Tuple_t &_tuple) const = 0;

	/** Writes a single table \a _table in the sqlite file.
	 *
	 * @param _table table to write to sqlite file
	 * @return true - everything ok, false - else
	 */
	virtual bool writeTable(const Table &_table) const = 0;

	/** Execute a given string SQL command \a _command.
	 *
	 * \param _command SQL statement to execute
	 * \return true - everything ok, false - command failed
	 */
	virtual bool executeSQLStatement(const std::string &_command) const = 0;

	/** Writes all tables in the database as an sqlite file.
	 *
	 * This takes the internal_database and writes it as a single
	 * table contained in a sqlite-compatible file.
	 *
	 * @return true - succesfully written, false - something went wrong
	 */
	virtual bool writeAllTables() const = 0;

	/** Getter for whether a database file was given.
	 *
	 * This is important for whether tables can be written to file and thus
	 * whether a unique rowid id (for parameter sets) can be obtained or not.
	 *
	 * @return true - filename given, false - else
	 */
	virtual bool isDatabaseFileGiven() const = 0;
};


#endif /* DATABASE_HPP_ */
