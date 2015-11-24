/*
 * TableDataOperator.hpp
 *
 *  Created on: Nov 24, 2015
 *      Author: heber
 */

#ifndef TABLEDATAOPERATOR_HPP_
#define TABLEDATAOPERATOR_HPP_

#include "BassoConfig.h"

#include "Database/Table.hpp"
#include "Database/types.hpp"

class Database;

/** This class represents the common stuff to all operations on table data
 * that pass along or accumulate information that needs extraction from
 * the table before.
 */
class TableDataOperator
{
public:
	//!> expose Table-internal to public part
	typedef std::vector<std::string> keys_t;

	/** Cstor of class InnerProblemDatabase.
	 *
	 * @param _db ref to database for obtaining tables
	 * @param _tablename name of table
	 * @param _accumulated_keys keys to accumulate
	 */
	TableDataOperator(
			const Database&_db);

	/** Cstor of class InnerProblemDatabase.
	 *
	 * @param _db ref to database for obtaining tables
	 * @param _tablename name of table
	 * @param _accumulated_keys_begin begin iteration for keys to accumulate
	 * @param _accumulated_keys_end end iteration for keys to accumulate
	 */
	TableDataOperator(
			const Database&_db,
			const std::string &_tablename,
			const keys_t::const_iterator &_accumulated_keys_begin,
			const keys_t::const_iterator &_accumulated_keys_end
			);

	/** Adds the accumulated information to the (hopefully) prepared \a _table.
	 *
	 * @param _table table to add accumulated information to
	 */
	void insertValues(
			Table &_table) const;

	typedef std::vector<
				std::pair<Database_types::types_t,
				Table::any_values_t> > valuevectors_t;

protected:
	// here we extend the TableDirectoryDatabase functionality

	/** Internal function to extract the tuples values into consistent
	 * vectors of values.
	 *
	 * @return vectors, one per \a accumulated_keys, each containing
	 * 			all the values over the tuples in the table.
	 */
	valuevectors_t extractValues() const;

	//!> accumulated data for placement in another table.
	mutable Table::TokenTypeMap_t accumulatedData;

	/** Internal function to accumulate the data into an internal map
	 * from the Table's Tuple_t's.
	 *
	 */
	virtual void extractData() const = 0;

	/** Internal function to average the data.
	 *
	 */
	virtual void finalizeData() const = 0;

	//!> reference to database to get tables from
	const Database& db;

	//!> name of table we accumulate data from
	const std::string tablename;

	//!> list of keys to accumulate
	Table::keys_t accumulated_keys;
};



#endif /* TABLEDATAOPERATOR_HPP_ */
