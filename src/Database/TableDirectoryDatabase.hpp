/*
 * TableDirectoryDatabase.hpp
 *
 *  Created on: Nov 24, 2015
 *      Author: heber
 */

#ifndef TABLEDIRECTORYDATABASE_HPP_
#define TABLEDIRECTORYDATABASE_HPP_

#include "BassoConfig.h"

#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "Database/TableDirectory.hpp"

/** This class implements a partial database where all tuples are stored in
 * a TableDirectory.
 */
class TableDirectoryDatabase : public Database
{
public:
	/** Cstor of class InnerProblemDatabase.
	 *
	 */
	TableDirectoryDatabase()
	{}

	/** Dstor of class InnerProblemDatabase.
	 *
	 */
	virtual ~TableDirectoryDatabase()
	{}

	/** Add a new table with name \a _name to this database.
	 *
	 * @param _name name of new table
	 * @return ref to the table, present one if it already existed
	 */
	Table& addTable(const std::string &_name)
	{ return directory.addTable(_name); }

	/** Removes a table with name \a _name from this database.
	 *
	 * @param _name name of new table
	 * @return true - table found and removed, false - else
	 */
	bool removeTable(const std::string &_name)
	{ return directory.removeTable(_name); }

	/** Removes contents of a table with name \a _name from this database.
	 *
	 * @param _name name of new table
	 * @return true - table found and contents removed, false - else
	 */
	virtual bool clearTable(const std::string &_name)
	{ return directory.clearTable(_name); }

	/** Getter for a table by its name \a _name.
	 *
	 * @param _name name of table
	 * @return ref to the table
	 */
	Table& getTable(const std::string &_name)
	{ return directory.getTable(_name); }

	/** Getter for a table by its name \a _name.
	 *
	 * @param _name name of table
	 * @return const ref to the table
	 */
	const Table& getTableConst(const std::string &_name) const
	{ return directory.getTableConst(_name); }

	/** Returns the number of tables present in the database.
	 *
	 * @return number of tables
	 */
	const size_t size() const
	{ return directory.size(); }

	/** Clears the whole database.
	 *
	 */
	void clear()
	{ directory.clear(); }

protected:
	//!> internal directory for all tables
	TableDirectory directory;
};



#endif /* TABLEDIRECTORYDATABASE_HPP_ */
