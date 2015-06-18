/*
 * Database_mock.hpp
 *
 *  Created on: Apr 15, 2015
 *      Author: heber
 */

#ifndef DATABASE_MOCK_HPP_
#define DATABASE_MOCK_HPP_

#include "Database.hpp"
#include "Table_mock.hpp"

class Database_mock : public Database
{
public:
	Database_mock() :
		table_mock("dummy_table")
	{}

	/** Add a new table with name \a _name to this database.
	 *
	 * @param _name name of new table
	 * @return ref to the table, present one if it already existed
	 */
	Table& addTable(const std::string &_name)
	{ return table_mock; }

	/** Removes a table with name \a _name from this database.
	 *
	 * @param _name name of new table
	 * @return true - table found and removed, false - else
	 */
	bool removeTable(const std::string &_name)
	{ return true; }

	/** Getter for a table by its name \a _name.
	 *
	 * @param _name name of table
	 * @return ref to the table
	 */
	Table& getTable(const std::string &_name)
	{ return table_mock; }

	/** Getter for a table by its name \a _name.
	 *
	 * @param _name name of table
	 * @return const ref to the table
	 */
	const Table& getTableConst(const std::string &_name) const
	{ return table_mock; }

	/** Returns the number of tables present in the database.
	 *
	 * @return number of tables
	 */
	const size_t size() const
	{ return 1; }

	/** Checks whether a tuple is present in the sqlite table without
	 * causing any warning if table is not uptodate.
	 *
	 * @param _table table to inspect
	 * @param _tuple tuple to look for
	 * @return true - tuple present, false - not present
	 */
	bool isTuplePresentInTable(
			const Table &_table,
			const Table::Tuple_t &_tuple) const
	{ return true; }

	/** Returns the rowid (sqlite specific primary key) for the requested
	 * \a _tuple in table \a _table.
	 *
	 * @param _table table to inspect
	 * @param _tuple tuple to look for
	 * @return rowid or (size_t)-1 of not present
	 */
	size_t getIdOfTuplePresentInTable(
			const Table &_table,
			const Table::Tuple_t &_tuple) const
	{ return 1; }

	/** Writes a single table \a _table in the sqlite file.
	 *
	 * @param _table table to write to sqlite file
	 * @return true - everything ok, false - else
	 */
	bool writeTable(const Table &_table) const
	{ return true; }

	/** Execute a given string SQL command \a _command.
	 *
	 * \param _command SQL statement to execute
	 * \return true - everything ok, false - command failed
	 */
	bool executeSQLStatement(const std::string &_command) const
	{ return true; }

	/** Writes all tables in the database as an sqlite file.
	 *
	 * This takes the internal_database and writes it as a single
	 * table contained in a sqlite-compatible file.
	 *
	 * @return true - succesfully written, false - something went wrong
	 */
	bool writeAllTables() const
	{ return true; }

private:
	//!> single mock table to return
	Table_mock table_mock;

};


#endif /* DATABASE_MOCK_HPP_ */
