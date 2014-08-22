/*
 * Table.hpp
 *
 *  Created on: Aug 22, 2014
 *      Author: heber
 */

#ifndef TABLE_HPP_
#define TABLE_HPP_

#include "BassoConfig.h"

#include "Database.hpp"

#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

class Database;

/** This class represents a single table in a Database.
 *
 */
class Table
{
public:
	/** Constructor for class Table.
	 *
	 */
	Table(const std::string &_tablename) :
		tablename(_tablename)
	{}

	/** Destructor for class Table.
	 *
	 */
	~Table() {}

	/** The Tuple_t contains the data as a map of key/value pairs.
	 */
	struct Tuple_t :
		public std::map<
			std::string,
			Database::typevariant_t >
	{
		bool operator<(const Tuple_t &_a) const;

		void replace(const std::string &_key, const Database::typevariant_t &_value);
	};

	/** Adds a data tuple to the table.
	 *
	 * @param _tuple tuple to add
	 */
	void addTuple(const Tuple_t &_tuple);

	/** Returns the number of tuples present in the table.
	 *
	 * @return number of tuples
	 */
	const size_t size() const
	{ return internal_table.size(); }

	/** Getter for the name of this Table.
	 *
	 * @return name of the table
	 */
	const std::string & getName() const
	{ return tablename; }

private:
	//!> grant Database access to private functions
	friend class Database;

	typedef std::set<std::string> keys_t;

	/** Returns the set of unique keys currently found in the table.
	 *
	 * @return set of unique keys
	 */
	keys_t getSetofUniqueKeys() const;

	typedef std::map< std::string, enum Database::types_t > KeyType_t;

	/** Determines for every key the unique type from the tuples present in
	 * the internal_table.
	 *
	 * @param _keys set of unique keys
	 * @return map from keys to types
	 */
	KeyType_t getKeyToTypeMap(
			const keys_t &_keys) const;

	/** Perform a sanity check whether each key has a value of the always
	 * same type over all tuples in internal_table.
	 *
	 * @param _keys set of unique keys to check
	 * @return true - all keys are ok, false - some keys have ambigious type
	 */
	bool checkTableSanity(const keys_t &_keys) const;

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
		for (internal_table_t::const_iterator iter = internal_table.begin();
				iter != internal_table.end(); ++iter) {
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

	//!> typedef for the format how the tuples are internally stored.
	typedef std::set<Tuple_t> internal_table_t;

	//!> internal table containing all tuples.
	internal_table_t internal_table;

	//!> contains the name of this table
	const std::string tablename;
};



#endif /* TABLE_HPP_ */
