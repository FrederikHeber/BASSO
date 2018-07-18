/*
 * Tuple_mock.hpp
 *
 *  Created on: Apr 15, 2015
 *      Author: heber
 */

#ifndef TUPLE_MOCK_HPP_
#define TUPLE_MOCK_HPP_

#include "Table.hpp"

class Table_mock;

/** This is a mock implementation of the tuple class that generates as little
 * overhead as possible.
 */
class Tuple_mock : public Table::Tuple_t
{
private:
	//!> grant Table_mock access to cstor and dstor
	friend class Table_mock;
	//!> grant Table_access to cstor and dstor
	friend class Table;

	Tuple_mock() {}

public:

	bool operator<(const Tuple_t &_a) const
	{ return true; }

	virtual ~Tuple_mock() {}

	/** Replaces a present (name, value) pair with a new \a _value.
	 *
	 * Asserts that name is already present.
	 *
	 * @param _key name of pair (column name)
	 * @param _value value of pair
	 */
	void replace(
			const std::string &_key,
			const Database_types::typevariant_t &_value)
	{}

	/** Inserts a new (name, value) pair.
	 *
	 * Asserts that name is not already present.
	 *
	 * @param _pair pair of (name,value)
	 * @param _type type of this pair
	 */
	void insert(
			const std::pair<
					std::string,
					Database_types::typevariant_t> &_pair,
			const enum Table::ColumnType _type)
	{}

	bool isParameter(const std::string &_name) const
	{ return false; }
};


#endif /* TUPLE_MOCK_HPP_ */
