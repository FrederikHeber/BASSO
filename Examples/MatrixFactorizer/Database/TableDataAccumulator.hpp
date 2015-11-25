/*
 * TableDataAccumulator.hpp
 *
 *  Created on: Nov 24, 2015
 *      Author: heber
 */

#ifndef TABLEDATAACCUMULATOR_HPP_
#define TABLEDATAACCUMULATOR_HPP_

#include "BassoConfig.h"

#include "MatrixFactorizer/Database/TableDataOperator.hpp"

class TableDataAccumulator : public TableDataOperator
{
public:
	/** Cstor of class TableDataAccumulator.
	 *
	 * @param _db ref to database for obtaining tables
	 */
	TableDataAccumulator(
			const Database&_db) :
		TableDataOperator(_db)
	{}

	/** Cstor of class TableDataAccumulator.
	 *
	 * @param _db ref to database for obtaining tables
	 * @param _tablename name of table
	 * @param _accumulated_keys_begin begin iteration for keys to accumulate
	 * @param _accumulated_keys_end end iteration for keys to accumulate
	 */
	TableDataAccumulator(
			const Database&_db,
			const std::string &_tablename,
			const keys_t::const_iterator &_accumulated_keys_begin,
			const keys_t::const_iterator &_accumulated_keys_end) :
		TableDataOperator(_db, _tablename, _accumulated_keys_begin, _accumulated_keys_end)
	{}

	/** Adds the required columns for the accumulated values to the \a _table.
	 *
	 * @param _table table to insert columns
	 * @param _type type of column
	 * @param _keyname key to insert (for each key four columns:
	 * 			min_, max_, avg_, and var_)
	 */
	static void prepareTableForAccumulatedValues(
			Table &_table,
			const Database_types::types_t &_type,
			const std::string &_keyname);

private:
	/** Internal function to accumulate the data into an internal map
	 * from the Table's Tuple_t's.
	 *
	 */
	void operateOnData() const;
};


#endif /* UNITTESTS_TABLEDATAACCUMULATOR_HPP_ */
