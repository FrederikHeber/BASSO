/*
 * Table.hpp
 *
 *  Created on: Aug 22, 2014
 *      Author: heber
 */

#ifndef TABLE_HPP_
#define TABLE_HPP_

#include "BassoConfig.h"

#include <boost/shared_ptr.hpp>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "Database/DefaultValue.hpp"
#include "Database/types.hpp"

class Database;

/** This class represents a single table in a Database.
 *
 */
class Table
{
protected:
	//!> grant Database access to private functions
	friend class Database;

	//!> typedef for this instance wrapped in a shared_ptr
	typedef boost::shared_ptr<Table> ptr;

	/** Constructor for class Table.
	 *
	 */
	Table(const std::string &_tablename) :
		tablename(_tablename),
		uptodate(true)
	{}

public:

	//!> enumeration on what a column represents
	enum ColumnType {
		Parameter=0,
		Data=1,
		MAX_ColumnType
	};

	//!> typedef for token of tuple and associated type
	typedef std::map<
			std::string,
			Database_types::typevariant_t > TokenTypeMap_t;

	/** The Tuple_t contains the data as a map of key/value pairs.
	 */
		class Tuple_t :
		public TokenTypeMap_t
	{
	protected:
			//!> grant only Table access to Tuple's cstor and dstor
			friend class Table;

			Tuple_t() {}
	public:

		virtual bool operator<(const Tuple_t &_a) const;

		/** Replaces a present (name, value) pair with a new \a _value.
		 *
		 * Asserts that name is already present.
		 *
		 * @param _key name of pair (column name)
		 * @param _value value of pair
		 */
		virtual void replace(
				const std::string &_key,
				const Database_types::typevariant_t &_value);

		/** Inserts a new (name, value) pair.
		 *
		 * Asserts that name is not already present.
		 *
		 * @param _pair pair of (name,value)
		 * @param _type type of this pair
		 */
		virtual void insert(
				const std::pair<
						std::string,
						Database_types::typevariant_t> &_pair,
				const enum ColumnType _type);

		virtual bool isParameter(const std::string &_name) const;

		/** Clears the internal map and the TypeMap.
		 *
		 */
		void clear();

	private:
		void insert(
				const std::pair<
						std::string,
						Database_types::typevariant_t> &_pair);

	private:
		typedef std::map<std::string, enum ColumnType> TypeMap_t;
		TypeMap_t TypeMap;
	};

	/** Getter for the default tuple associated to this table.
	 *
	 * @return ref to tuple
	 */
	virtual Tuple_t& getTuple()
	{ return internal_tuple; }

	/** Adds a data tuple to the table.
	 *
	 * @param _tuple tuple to add
	 */
	virtual void addTuple(const Tuple_t &_tuple);

	/** Returns the number of tuples present in the table.
	 *
	 * @return number of tuples
	 */
	virtual const size_t size() const
	{ return internal_table.size(); }

	/** Getter for the name of this Table.
	 *
	 * @return name of the table
	 */
	virtual const std::string & getName() const
	{ return tablename; }

	/** Getter whether table has not been changed since last update.
	 *
	 * @return true - table is unchanged, false - else
	 */
	virtual const bool isUptodate() const
	{ return uptodate; }

	/** Removes all tuples contained in a table.
	 *
	 */
	virtual void clear();

	/** Returns whether table is empty.
	 *
	 * @return true - no tuples present, false - else
	 */
	virtual bool empty() const
	{ return internal_table.empty(); }

private:

	//!> unique set of all column names
	typedef std::set<std::string> keys_t;

	/** Returns the set of unique keys currently found in the table.
	 *
	 * @return set of unique keys
	 */
	keys_t getSetofUniqueKeys() const;

	//!> map of keys (column names) to types (all possible type variants)
	typedef std::map< std::string, enum Database_types::types_t > KeyType_t;

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
	 * @return true - all keys are ok, false - some keys have ambiguous type
	 */
	bool checkTableSanity(const keys_t &_keys) const;

	//!> typedef for a single row in stringized form (vector of strings)
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


	/** Returns all tuples as an array (vector of vectors) of string
	 * values.
	 *
	 * \param _KeyTypes map with type for each key (column name)
	 * \param _allowed_keys vector of all allowed keys (filter for column names)
	 * \return vector of tuples, each a vector of string values of each cell
	 * 		   in the table
	 */
	std::vector<values_t> convertTuplesToValueVector(
			const KeyType_t &_KeyTypes,
			const keys_t &_allowed_keys) const;

private:

	//!> typedef for the format how the tuples are internally stored.
	typedef std::set<Tuple_t> internal_table_t;

	//!> internal table containing all tuples.
	internal_table_t internal_table;

	//!> contains the name of this table
	const std::string tablename;

	/** Flag indicating whether written SQL database is up-to-date
	 * with respect to this table.
	 */
	mutable bool uptodate;

	/** Internal tuple used for this table.
	 *
	 */
	Tuple_t internal_tuple;
};



#endif /* TABLE_HPP_ */
