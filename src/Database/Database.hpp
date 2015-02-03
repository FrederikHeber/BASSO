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
#include <boost/variant.hpp>
#include <map>
#include <string>
#include <vector>

#include "Database/DefaultValue.hpp"

class Table;

/** The Database class provides a single instance to contain all iteration-
 * related data. The data is eventually placed in a single-file (sqlite)
 * allowing for command-line parsing and data extraction for creating figures.
 */
class Database
{
public:
	//!> typedef for wrapping Database instance in shared_ptr
	typedef boost::shared_ptr<Database> Database_ptr_t;

	typedef boost::variant<int, double, std::string > typevariant_t;

	enum types_t
	{
		inttype=0,
		doubletype=1,
		valchartype=2,
		MAX_TYPES
	};

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
	bool writeSQLitefile();

private:
	//!> states whether a database file name has been given or not.
	bool DatabaseFileGiven;

	//!> states whether we replace already present parameter tuples or just add new ones
	bool ReplacePresentParameterTuples;

	//!> contains the filename of the database for storing
	std::string filename;

	//!> static vector with all type names
	static std::vector<std::string> TypeNames;

	//!> typedef for the set 1of tables
	typedef std::map<std::string, Table> tables_t;

	//!> list of all tables of this database
	tables_t tables;

	//!> maximum number of keys per table
	const unsigned int MaxKeys;
};


#endif /* DATABASE_HPP_ */
