/*
 * Table_mock.hpp
 *
 *  Created on: Apr 15, 2015
 *      Author: heber
 */

#ifndef TABLE_MOCK_HPP_
#define TABLE_MOCK_HPP_

#include "Table.hpp"
#include "Tuple_mock.hpp"

/** This is a mock implementation of the Table class that generates as little
 * overhead as possible.
 */
class Table_mock : public Table
{
public:
	Table_mock(const std::string &_name) :
		Table(_name)
	{}

	/** Getter for the default tuple associated to this table.
	 *
	 * @return ref to tuple
	 */
	Tuple_t& getTuple()
	{ return dummy_tuple; }

	/** Returns the number of tuples present in the table.
	 *
	 * @return number of tuples
	 */
	const size_t size() const
	{ return 1; }

	/** Getter whether table has not been changed since last update.
	 *
	 * @return true - table is unchanged, false - else
	 */
	const bool isUptodate() const
	{ return true; }

	/** Removes all tuples contained in a table.
	 *
	 */
	void clear()
	{ }

	/** Returns whether table is empty.
	 *
	 * @return true - no tuples present, false - else
	 */
	bool empty() const
	{ return false; }

private:
	//!> mock tuple that does nothing
	Tuple_mock dummy_tuple;
};


#endif /* TABLE_MOCK_HPP_ */
