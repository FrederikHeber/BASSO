/*
 * TableDirectory.hpp
 *
 *  Created on: Nov 23, 2015
 *      Author: heber
 */

#ifndef TABLEDIRECTORY_HPP_
#define TABLEDIRECTORY_HPP_

#include "BassoConfig.h"

#include <map>

#include "Database/Table.hpp"

class Database;

/** This structure maintains the set of stables as a directory,
 * i.e. each table is associated with a name.
 */
struct TableDirectory
{
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

	/** Clears the whole directory.
	 *
	 */
	void clear()
	{ tables.clear(); }

private:
	//!> typedef for the set of tables
	typedef std::map<std::string, Table::ptr> tables_t;

	//!> list of all tables of this database
	tables_t tables;

public:
	//!> typedef to give const iterator to outside
	typedef tables_t::const_iterator TableIterator_t;

	/** Get begin iterator to tables.
	 *
	 * @return begin iterator
	 */
	TableIterator_t begin() const
	{ return tables.begin(); }

	/** Get end iterator to tables.
	 *
	 * @return end iterator
	 */
	TableIterator_t end() const
	{ return tables.end(); }
};



#endif /* TABLEDIRECTORY_HPP_ */
