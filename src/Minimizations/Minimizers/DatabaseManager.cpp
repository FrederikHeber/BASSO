/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*
 * GeneralMinimizer_DatabaseManager.cpp
 *
 *  Created on: Nov 25, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "GeneralMinimizer.hpp"

#include <cassert>
#include <cmath>

#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

#include "Database/Database.hpp"
#include "Log/Logging.hpp"

static void NoOpAddParams(Table::Tuple_t &, const bool)
{}


DatabaseManager::DatabaseManager(Database &_database) :
	database(_database),
	add_param_callback(boost::bind(NoOpAddParams, _1, _2)),
	parameter_key(0),
	parameters_table(database.addTable("parameters")),
	data_per_iteration_table(database.addTable("data_per_iteration")),
	data_overall_table(database.addTable("data_overall"))
{}

size_t DatabaseManager::size() const
{
	return database.size();
}

void DatabaseManager::setParameterKey(
		double _val_NormX,
		double _val_NormY,
		const unsigned int _N,
		const unsigned int _dim,
		const int _MaxOuterIterations) const
{
	// convert tuple values for "inf" case (not valid SQL number)
	if (isinf(_val_NormX))
		_val_NormX = 0.;
	if (isinf(_val_NormY))
		_val_NormY = 0.;

	// create tuple
	Table::Tuple_t &parameter_tuple = parameters_table.getTuple();
	const bool do_replace = !parameter_tuple.empty();
	if (!do_replace) {
		parameter_tuple.insert( std::make_pair("p", _val_NormX), Table::Parameter);
		parameter_tuple.insert( std::make_pair("r", _val_NormY), Table::Parameter);
		parameter_tuple.insert( std::make_pair("N", (int)_N), Table::Parameter);
		parameter_tuple.insert( std::make_pair("dim", (int)_dim), Table::Parameter);
		parameter_tuple.insert( std::make_pair("max_iterations", _MaxOuterIterations), Table::Parameter);
		// add additional parameters from derived classes
		add_param_callback(parameter_tuple, do_replace);
		// add additional parameters specified from user
		for (std::vector<std::string>::const_iterator iter = tuple_params.begin();
				iter != tuple_params.end(); ) {
			const std::string &token = (*iter++);
			const std::string &value = (*iter++);
			try {
				const int int_value = boost::lexical_cast<int>(value);
				parameter_tuple.insert( std::make_pair(token, int_value), Table::Parameter);
			} catch(const boost::bad_lexical_cast &) {
				try {
					const double double_value = boost::lexical_cast<double>(value);
					parameter_tuple.insert( std::make_pair(token, double_value), Table::Parameter);
				} catch(const boost::bad_lexical_cast &) {
					parameter_tuple.insert( std::make_pair(token, value), Table::Parameter);
				}
			}
		}
	} else {
		parameter_tuple.replace( "p", _val_NormX);
		parameter_tuple.replace( "r", _val_NormY);
		parameter_tuple.replace( "N", (int)_N);
		parameter_tuple.replace( "dim", (int)_dim);
		parameter_tuple.replace( "max_iterations", _MaxOuterIterations);
		// add additional parameters from derived classes
		add_param_callback(parameter_tuple, do_replace);
		// add additional parameters specified from user
		for (std::vector<std::string>::const_iterator iter = tuple_params.begin();
				iter != tuple_params.end(); ) {
			const std::string &token = (*iter++);
			const std::string &value = (*iter++);
			try {
				const int int_value = boost::lexical_cast<int>(value);
				const double double_value = boost::lexical_cast<double>(value);
				if ((double)int_value == double_value)
					parameter_tuple.replace( token, int_value);
				else
					parameter_tuple.replace( token, double_value);
			} catch(const boost::bad_lexical_cast &) {
				parameter_tuple.replace( token, value);
			}
		}
	}
	// we need to add it, otherwise we cannot use checks for presence
	// as they rely on complete table info
	parameters_table.addTuple(parameter_tuple);

	size_t rowid = 0;
	if (database.isDatabaseFileGiven()) {
		// check for presence
		if (!database.isTuplePresentInTable(parameters_table, parameter_tuple)) {
			LOG(debug, "Parameter tuple not present, adding to table.");
			database.writeTable(parameters_table);
		}
		// and store rowid
		rowid = database.getIdOfTuplePresentInTable(
					parameters_table, parameter_tuple);
		// clear table such that present tuple is not stored again
		database.clearTable(parameters_table.getName());
		LOG(debug, "Setting parameter_key to " << rowid);
	} else {
		// else set rowid to arbitrary value as there is no file anyway
		rowid = 1;
	}
	const_cast<size_t &>(parameter_key) = rowid;
}

Table::Tuple_t & DatabaseManager::preparePerIterationTuple() const
{
	assert(parameter_key != 0);

	Table::Tuple_t &per_iteration_tuple = data_per_iteration_table.getTuple();
	per_iteration_tuple.insert( std::make_pair("parameters_fk", (int)parameter_key), Table::Parameter);
	per_iteration_tuple.insert( std::make_pair("iteration", (int)0), Table::Data);
	per_iteration_tuple.insert( std::make_pair("stepwidth", 0.), Table::Data);
	per_iteration_tuple.insert( std::make_pair("residual", 0.), Table::Data);
	per_iteration_tuple.insert( std::make_pair("relative_residual", 0.), Table::Data);
	per_iteration_tuple.insert( std::make_pair("error", 0.), Table::Data);
	per_iteration_tuple.insert( std::make_pair("bregman_distance", 0.), Table::Data);
	per_iteration_tuple.insert( std::make_pair("updated_index", (int)0), Table::Data);
	return per_iteration_tuple;
}

Table::Tuple_t & DatabaseManager::prepareOverallTuple() const
{
	assert(parameter_key != 0);

	Table::Tuple_t &overall_tuple = data_overall_table.getTuple();
	overall_tuple.insert( std::make_pair("parameters_fk", (int)parameter_key), Table::Parameter);
	overall_tuple.insert( std::make_pair("iterations", (int)0), Table::Data);
	overall_tuple.insert( std::make_pair("residual", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("relative_residual", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("runtime", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("element_creation_operations", (int)0), Table::Data);
	overall_tuple.insert( std::make_pair("linear_time_operations", (int)0), Table::Data);
	overall_tuple.insert( std::make_pair("quadratic_time_operations", (int)0), Table::Data);
	overall_tuple.insert( std::make_pair("element_creation_runtime", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("linear_time_runtime", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("quadratic_time_runtime", 0.), Table::Data);
	return overall_tuple;
}

void DatabaseManager::finalizeOverallTuple(
		Table::Tuple_t &_overall_tuple,
		QuickAccessReferences &_refs) const
{
	_overall_tuple.replace( "element_creation_operations",
			(int)(_refs.SpaceX.getOpCounts().getTotalConstantCounts()
					+_refs.SpaceY.getOpCounts().getTotalConstantCounts()
					+_refs.DualSpaceX.getOpCounts().getTotalConstantCounts()
					+_refs.DualSpaceY.getOpCounts().getTotalConstantCounts()));
	_overall_tuple.replace( "linear_time_operations",
			(int)(_refs.SpaceX.getOpCounts().getTotalLinearCounts()
					+_refs.SpaceY.getOpCounts().getTotalLinearCounts()
					+_refs.DualSpaceX.getOpCounts().getTotalLinearCounts()
					+_refs.DualSpaceY.getOpCounts().getTotalLinearCounts()));
	_overall_tuple.replace( "quadratic_time_operations",
			(int)(_refs.A.getCount()+_refs.A_t.getCount()) );
	// NOTE: due to Eigen's lazy evaluation runtime is not measured accurately
	_overall_tuple.replace( "element_creation_runtime",
			boost::chrono::duration<double>(
					_refs.SpaceX.getOpCounts().getTotalConstantTimings()
					+_refs.SpaceY.getOpCounts().getTotalConstantTimings()
					+_refs.DualSpaceX.getOpCounts().getTotalConstantTimings()
					+_refs.DualSpaceY.getOpCounts().getTotalConstantTimings()).count());
	_overall_tuple.replace( "linear_time_runtime",
			boost::chrono::duration<double>(
					_refs.SpaceX.getOpCounts().getTotalLinearTimings()
					+_refs.SpaceY.getOpCounts().getTotalLinearTimings()
					+_refs.DualSpaceX.getOpCounts().getTotalLinearTimings()
					+_refs.DualSpaceY.getOpCounts().getTotalLinearTimings()).count());
	_overall_tuple.replace( "quadratic_time_runtime",
			boost::chrono::duration<double>(
					_refs.A.getTiming()+_refs.A_t.getTiming()).count() );
	// NOTE: due to Eigen's lazy evaluation runtime is not measured accurately
}

bool DatabaseManager::createViews() const
{
	// write tables beforehand
	database.writeAllTables();

	// we create views overall, and per_iteration (which were present
	// before the switch to distinct data and parameters table)
	bool status = true;
	{
		// check whether tables are present and contain elements
		status &= !data_overall_table.empty();
		status &= !data_per_iteration_table.empty();
		// we don't  check for parameter table's non-emptiness as is
		// possibly might be if the used parameter tuple is already in
		// the database, see setParameterKey()
	}
	if (!status) {
		LOG(error, "(Some of the) Required Tables are empty, not creating views.");
	}
	if (status) {
		std::stringstream sql;
		sql << "CREATE VIEW IF NOT EXISTS overall AS SELECT * FROM parameters p INNER JOIN data_overall d ON p.rowid = d.parameters_fk";
		LOG(trace, "SQL: " << sql.str());
		status &= database.executeSQLStatement(sql.str());
	}
	if (status) {
		std::stringstream sql;
		sql << "CREATE VIEW IF NOT EXISTS per_iteration AS SELECT * FROM parameters p INNER JOIN data_per_iteration d ON p.rowid = d.parameters_fk";
		LOG(trace, "SQL: " << sql.str());
		status &= database.executeSQLStatement(sql.str());
	}
	return status;
}

void DatabaseManager::setAdditionalTupleParameters(
			const std::vector<std::string> &_tuple_params)
{
	assert ( _tuple_params.size() % 2 == 0);
	const_cast< std::vector<std::string> &>(tuple_params) = _tuple_params;
}
