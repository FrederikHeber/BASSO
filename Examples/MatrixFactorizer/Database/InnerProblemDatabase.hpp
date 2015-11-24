/*
 * InnerProblemDatabase.hpp
 *
 *  Created on: Nov 23, 2015
 *      Author: heber
 */

#ifndef INNERPROBLEMDATABASE_HPP_
#define INNERPROBLEMDATABASE_HPP_

#include "BassoConfig.h"

#include "Database/TableDirectoryDatabase.hpp"

class InnerProblemDatabase : public TableDirectoryDatabase
{
public:
	/** Cstor of class InnerProblemDatabase.
	 *
	 */
	InnerProblemDatabase();

	/** Dstor of class InnerProblemDatabase.
	 *
	 */
	virtual ~InnerProblemDatabase();

	void setDatabaseFile( const std::string &_filename)
	{}

	/** Setter for whether tuples are replaced (with respect to their parameter
	 * part or we just add new ones regardless of same-parametered present ones.
	 *
	 * \warning Replacement can take quite some time for a larger database
	 * files (e.g. >10Mb ~ secs, >100Mb ~ mins).
	 *
	 * @param _flag true - replace, false - just add
	 */
	void setReplacePresentParameterTuples(const bool _flag)
	{}

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

	/** Getter for whether a database file was given.
	 *
	 * This is important for whether tables can be written to file and thus
	 * whether a unique rowid id (for parameter sets) can be obtained or not.
	 *
	 * @return true - filename given, false - else
	 */
	bool isDatabaseFileGiven() const
	{ return true; }

private:

	/** Addes the required columns for the accumulated values to the \a _table.
	 *
	 * @param _table table to insert columns
	 * @param _keyname key to insert (for each key four columns:
	 * 			min_, max_, avg_, and var_)
	 */
	static void prepareTableForAccumulatedValues(
			Table &_table,
			const std::string &_keyname);

	/** Addes the accumulated information to the (hopefully) prepared \a _table.
	 *
	 * @param _table table to add accumulated information to
	 */
	void insertAccumulatedValues(
			Table &_table) const;

public:
	// here we extend the TableDirectoryDatabase functionality

	void accumulateData() const;

	//!> accumulated data for placement in another table.
	mutable Table::TokenTypeMap_t accumulatedData;
};


#endif /* INNERPROBLEMDATABASE_HPP_ */
