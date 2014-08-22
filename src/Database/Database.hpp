/*
 * Database.hpp
 *
 *  Created on: Jun 12, 2014
 *      Author: heber
 */

#ifndef DATABASE_HPP_
#define DATABASE_HPP_

#include "BassoConfig.h"

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

	//!> contains the filename of the database for storing
	std::string filename;

	//!> static vector with all type names
	static std::vector<std::string> TypeNames;

	//!> typedef for the set 1of tables
	typedef std::map<std::string, Table> tables_t;

	//!> list of all tables of this database
	tables_t tables;

	//!> maximum number of keys per table
	const unsigned int MAXKEYS;
};


#endif /* DATABASE_HPP_ */
