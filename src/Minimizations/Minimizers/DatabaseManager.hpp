/*
 * DatabaseManager.hpp
 *
 *  Created on: Nov 25, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_DATABASEMANAGER_HPP_
#define MINIMIZATIONS_MINIMIZERS_DATABASEMANAGER_HPP_

#include "BassoConfig.h"

#include <string>
#include <vector>

#include <boost/function.hpp>

#include "Database/Table.hpp"
#include "Minimizations/InverseProblems/QuickAccessReferences.hpp"

class Database;

struct DatabaseManager
{
	//!> typedef for callback function
	typedef boost::function<
			void (Table::Tuple_t &, const bool)> addAdditionalParameters_t;

	/** Cstor of struct DatabaseManager.
	 *
	 * @param _database database to manage
	 */
	DatabaseManager(Database &_database);

	/** Returns the number of tables in the database.
	 *
	 * @return number of tables
	 */
	size_t size() const;

	/** Create (deprecated) overall and per_iteration tables as views.
	 *
	 * \return true - views created, false - statements failed
	 */
	bool createViews() const;

	/** Setter for the callback function to add more parameter columns.
	 *
	 * @param _add_param_callback call back to set
	 */
	void setAddParamsCallback(
			const addAdditionalParameters_t _add_param_callback)
	{ const_cast<addAdditionalParameters_t&>(add_param_callback) =
			_add_param_callback;	}

	/** Sets the vector with additional parameters to make tuples unique in database.
	 *
	 * @param _tuple_params vector of strings in pairs of two ("key" and "value")
	 */
	void setAdditionalTupleParameters(
			const std::vector<std::string> &_tuple_params);

	/** Sets the parameter key by supplying a given tuple of values
	 *
	 * @param _val_NormX norm of space X
	 * @param _val_NormY norm of space Y
	 * @param _N number of search directions
	 * @param _dim dimension of solution vector
	 * @param _MaxOuterIterations maximum number of outer iterations
	 */
	void setParameterKey(
			double _val_NormX,
			double _val_NormY,
			const unsigned int _N,
			const unsigned int _dim,
			const int _MaxOuterIterations) const;

	Table::Tuple_t & preparePerIterationTuple() const;

	Table::Tuple_t & prepareOverallTuple() const;

	void finalizeOverallTuple(
			Table::Tuple_t &_overall_tuple,
			QuickAccessReferences &_refs) const;

	/** reference to an external database where we store infomation
	 * about the behavior of the iteration procedure.
	 */
	Database &database;

	//!> callback to use when parameters table is created
	const addAdditionalParameters_t add_param_callback;

	//!> additional parameter, value pairs that are added to each submitted tuple
	const std::vector<std::string> tuple_params;
	//!> contains the primary key of the current parameter set
	const size_t parameter_key;

	//!> table with parameter information
	Table& parameters_table;
	//!> table with per iteration step information
	Table& data_per_iteration_table;
	//!> table with overall iteration information
	Table& data_overall_table;

};

#endif /* MINIMIZATIONS_MINIMIZERS_DATABASEMANAGER_HPP_ */
