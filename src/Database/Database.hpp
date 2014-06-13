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
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "Database/DefaultValue.hpp"

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

	/** The Tuple_t contains the data as a map of key/value pairs.
	 */
	struct Tuple_t :
		public std::map<
			std::string,
			typevariant_t >
	{
		bool operator<(const Tuple_t &_a) const;

		void replace(const std::string &_key, const typevariant_t &_value);
	};

	/** Adds a data tuple to the database.
	 *
	 * @param _tuple tuple to add
	 */
	void addTuple(const Tuple_t &_tuple);

	/** Returns the number of tuples present in the database.
	 *
	 * @return number of tuples
	 */
	const size_t size() const
	{ return internal_database.size(); }

	void setDatabaseFile( const std::string &_filename)
	{ filename = _filename; DatabaseFileGiven = true; }

private:
	/** Internal function to write the database as an sqlite file.
	 *
	 * This takes the internal_database and writes it as a single
	 * table contained in a sqlite-compatible file.
	 *
	 * @return true - succesfully written, false - something went wrong
	 */
	bool writeSQLitefile();

	typedef std::set<std::string> keys_t;

	/** Returns the set of unique keys currently found in the database.
	 *
	 * @return set of unique keys
	 */
	keys_t getSetofUniqueKeys() const;

	typedef std::map< std::string, enum types_t > KeyType_t;

	/** Determines for every key the unique type from the tuples present in
	 * the internal_database.
	 *
	 * @param _keys set of unique keys
	 * @return map from keys to types
	 */
	KeyType_t getKeyToTypeMap(
			const keys_t &_keys) const;

	/** Perform a sanity check whether each key has a value of the always
	 * same type over all tuples in internal_database.
	 *
	 * @param _keys set of unique keys to check
	 * @return true - all keys are ok, false - some keys have ambigious type
	 */
	bool checkDatabaseSanity(const keys_t &_keys) const;

    typedef std::vector< std::string > values_t;

    /** Returns all values of \a _keys for each tuple.
     *
     * @param _key key to look for
     * @return vector of values of the corresponding keys
     */
    template <class T>
    values_t getAllValuesPerType(const std::string &_key) const
    {
    	values_t values;
		for (internal_database_t::const_iterator iter = internal_database.begin();
				iter != internal_database.end(); ++iter) {
			Tuple_t::const_iterator finditer = iter->find(_key);
			// if key present in given ones
			DefaultValue<T> defaultval;
			T value = defaultval.get();
			if (finditer != iter->end()) {
				value = boost::get<T>(finditer->second);
			}
			// append value
			{
				std::stringstream output;
				output << value;
				values.push_back(output.str());
			}
		}
    	return values;
    }

private:
	//!> states whether a database file name has been given or not.
	bool DatabaseFileGiven;

	//!> contains the filename of the database for storing
	std::string filename;

	//!> typedef for the format how the tuples are internally stored.
	typedef std::set<Tuple_t> internal_database_t;

	//!> internal database containing all tuples.
	internal_database_t internal_database;

	//!> static vector with all type names
	static std::vector<std::string> TypeNames;
};


#endif /* DATABASE_HPP_ */
