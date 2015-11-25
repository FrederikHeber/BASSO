/*
 * IterationInformation.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_DATABASE_ITERATIONINFORMATION_HPP_
#define MATRIXFACTORIZERBASE_DATABASE_ITERATIONINFORMATION_HPP_

#include "BassoConfig.h"

class MatrixFactorizerOptions;

#include "Log/Logging.hpp"

#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "Minimizations/types.hpp"

/** This stores the connect to a database and required tables
 *
 */
struct IterationInformation
{
	/** Constructor for class IterationInformation.
	 *
	 * Initializes database, tables, and tuples.
	 *
	 * @param _opts factorization options
	 * @param _innersize inner size of matrix to factorize
	 * @param _outersize outer size of matrix to factorize
	 */
	IterationInformation(
			const MatrixFactorizerOptions &_opts,
			const unsigned int _innersize,
			const unsigned int _outersize);

	/** Destructor for class IterationInformation.
	 *
	 * This creates views after (some) data has been written to.
	 *
	 */
	~IterationInformation();

	//!> enumeration of all tables contained in this struct
	enum TableType {
		LoopTable,
		OverallTable
	};

	/** Replaces a value of key \a _keyname in the tuple of the respective
	 * matrix \a _table.
	 *
	 * \param _type type of table to add to.
	 * \param _keyname name of key whose value to replace
	 * \param _value value to replace
	 */
	template <typename T>
	void replace(
			const enum TableType _table,
			const std::string &_keyname,
			const T _value
			)
	{
		switch(_table) {
		case LoopTable:
			loop_tuple.replace(_keyname, _value);
			break;
		case OverallTable:
			overall_tuple.replace(_keyname, _value);
			break;
		default:
			BOOST_LOG_TRIVIAL(error)
				<< "Unknown table in IterationInformation::replace()";
		}
	}

	/** Adds the currently present data tuple as it is for the respective
	 * matrix \a _table to the table.
	 *
	 * \param _type type of table to add to.
	 */
	void addTuple(
			const enum TableType _table
			);

	/** Prepares the parameters table in the database.
	 *
	 * @param _innerSize inner size of problem matrix
	 * @param _outerSize outer size of problem matrix
	 * @param _opts options with other information for setting parameters entry
	 * @return rowid of the set parameter key
	 */
	size_t prepareParametersTable(
			const size_t _innerSize,
			const size_t _outerSize,
			const MatrixFactorizerOptions &_opts
			);

	/** Getter to \a loop_table.
	 *
	 * @return ref to loop_table
	 */
	Table& getLoopTable()
	{ return loop_table; }

	/** Getter to \a loop_overall_table.
	 *
	 * @return ref to loop_overall_table
	 */
	Table& getOverallTable()
	{ return loop_overall_table; }

private:
	/** Creates views to maintain compatibility with earlier database design.
	 *
	 * @return true - success
	 */
	bool createViews();

	/** Prepares a dummy "overall" table used for the inner optimization
	 * values.
	 *
	 * @param _dummytable table to use
	 * @param _parameter_key parameter key
	 * @return tuple for this table
	 */
	static Table::Tuple_t & prepareOverallTuple(
			Table &_dummytable,
			const int _parameter_key);

private:
	//!> database connection
	Database_ptr_t database;
	//!> table with per loop information
	Table& loop_table;
	//!> table with overall minimization information
	Table& loop_overall_table;
	//!> data tuple for loop_table
	Table::Tuple_t& loop_tuple;
	//!> data tuple for loop_overall_table
	Table::Tuple_t& overall_tuple;
};



#endif /* MATRIXFACTORIZERBASE_DATABASE_ITERATIONINFORMATION_HPP_ */
