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
 * IterationInformation.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "IterationInformation.hpp"

#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

#include "Database/TableDirectoryDatabase_mock.hpp"
#include "Log/Logging.hpp"
#include "MatrixFactorizer/Database/TableDataAccumulator.hpp"
#include "MatrixFactorizer/Options/MatrixFactorizerOptions.hpp"
#include "Minimizations/Minimizers/DatabaseManager.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"

void addTupleParametersToTuple(
		const std::vector<std::string> &_tuple_parameters,
		Table::Tuple_t& _parameter_tuple
		)
{
	for (std::vector<std::string>::const_iterator iter = _tuple_parameters.begin();
			iter != _tuple_parameters.end(); ) {
		const std::string &token = (*iter++);
		const std::string &value = (*iter++);
		try {
			const int int_value = boost::lexical_cast<int>(value);
			_parameter_tuple.insert( std::make_pair(token, int_value), Table::Parameter);
			LOG(debug, " Adding additional integer parameter (" << token << "," << int_value << ") to loop tuple.");
		} catch(const boost::bad_lexical_cast &) {
			try {
				const double double_value = boost::lexical_cast<double>(value);
				_parameter_tuple.insert( std::make_pair(token, double_value), Table::Parameter);
				LOG(debug, " Adding additional double parameter (" << token << "," << double_value << ") to loop tuple.");
			} catch(const boost::bad_lexical_cast &) {
				_parameter_tuple.insert( std::make_pair(token, value), Table::Parameter);
				LOG(debug, " Adding additional string parameter (" << token << "," << value << ") to loop tuple.");
			}
		}
	}
}

IterationInformation::IterationInformation(
		const MatrixFactorizerOptions &_opts,
		const unsigned int _innersize,
		const unsigned int _outersize) :
			database(SolverFactory::createDatabase(_opts)),
			loop_table(database->addTable("data_loop")),
			loop_overall_table(database->addTable("data_loop_overall")),
			loop_tuple(loop_table.getTuple()),
			overall_tuple(loop_overall_table.getTuple())
{
	/// install parameter set and get its key
	const size_t parameter_key = prepareParametersTable(
			_innersize,
			_outersize,
			_opts);

	loop_tuple.insert( std::make_pair("parameters_fk", (int)parameter_key), Table::Parameter);
	loop_tuple.insert( std::make_pair("loop_nr", (int)0), Table::Data);
	loop_tuple.insert( std::make_pair("residual", 0.), Table::Data);
	loop_tuple.insert( std::make_pair("scaling_change", 0.), Table::Data);

	/// now added columns for accumulate values
	{
		// we need to get at the correct type for each desired accumulate
		// value. The type is contained in Minimizer's DatabaseManager
		// class. Hence, we instantiate it here, get the types and drop it.
		// in order to get each columns type, we instantiate a GeneralMinimizer
		Database_ptr_t tmpdatabase(new TableDirectoryDatabase_mock);
		DatabaseManager manager(*tmpdatabase);
		manager.setParameterKey(2.,2.,1,1,0);
		Table::Tuple_t &overall_tuple = manager.prepareOverallTuple();
		for (std::vector<std::string>::const_iterator iter = _opts.overall_keys.begin();
				iter != _opts.overall_keys.end(); ++iter) {
			// find key in overall_tuple
			Table::TokenTypeMap_t::const_iterator tupleiter =
					overall_tuple.find(*iter);
			if (tupleiter == overall_tuple.end()) {
				LOG(error, "Cannot find key " << *iter << "for accumulated values in overall table columns.");
			} else {
				TableDataAccumulator::prepareTableForAccumulatedValues(
						loop_table,
						Table::getVariantsType(tupleiter->second),
						(*iter)+"_projection");
				TableDataAccumulator::prepareTableForAccumulatedValues(
						loop_table,
						Table::getVariantsType(tupleiter->second),
						(*iter)+"_minimization");
			}
		}
	}

	overall_tuple.insert( std::make_pair("parameters_fk", (int)parameter_key), Table::Parameter);
	overall_tuple.insert( std::make_pair("loops", (int)0), Table::Data);
	overall_tuple.insert( std::make_pair("residual", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("runtime", 0.), Table::Data);
}

IterationInformation::~IterationInformation()
{
	// add views after tables have been created and filled
	if (!createViews()) {
		LOG(warning, "Could not create overall or per_iteration views.");
	}
}

void IterationInformation::addTuple(
		const enum TableType _table
		)
{
	switch(_table) {
	case LoopTable:
		loop_table.addTuple(loop_tuple);
		break;
	case OverallTable:
		loop_overall_table.addTuple(overall_tuple);
		break;
	default:
		LOG(error, "Unknown table in IterationInformation::addTuple()");
	}
}

Table::Tuple_t & IterationInformation::prepareOverallTuple(
		Table &_dummytable,
		const int _parameter_key
		)
{
	assert(_parameter_key != 0);

	Table::Tuple_t &dummy_tuple = _dummytable.getTuple();
	dummy_tuple.insert( std::make_pair("parameters_fk", (int)_parameter_key), Table::Parameter);
	dummy_tuple.insert( std::make_pair("iterations", (int)0), Table::Data);
	dummy_tuple.insert( std::make_pair("residual", 0.), Table::Data);
	dummy_tuple.insert( std::make_pair("relative_residual", 0.), Table::Data);
	dummy_tuple.insert( std::make_pair("runtime", 0.), Table::Data);
	dummy_tuple.insert( std::make_pair("element_creation_operations", (int)0), Table::Data);
	dummy_tuple.insert( std::make_pair("linear_time_operations", (int)0), Table::Data);
	dummy_tuple.insert( std::make_pair("quadratic_time_operations", (int)0), Table::Data);
	dummy_tuple.insert( std::make_pair("element_creation_runtime", 0.), Table::Data);
	dummy_tuple.insert( std::make_pair("linear_time_runtime", 0.), Table::Data);
	dummy_tuple.insert( std::make_pair("quadratic_time_runtime", 0.), Table::Data);
	return dummy_tuple;
}

bool IterationInformation::createViews()
{
	// write tables beforehand
	database->writeAllTables();

	bool status = true;
	{
		// check whether tables are present and contain elements
		status &= !loop_table.empty();
		status &= !loop_overall_table.empty();
		// we don't  check for parameter table's non-emptiness as is
		// possibly might be if the used parameter tuple is already in
		// the database, see setParameterKey()
	}
	if (!status) {
		LOG(error, "(Some of the) Required Tables are empty, not creating views.");
	}
	if (status) {
		std::stringstream sql;
		sql << "CREATE VIEW IF NOT EXISTS loop AS SELECT * FROM parameters p INNER JOIN data_loop d ON p.rowid = d.parameters_fk";
		LOG(trace, "SQL: " << sql.str());
		status &= database->executeSQLStatement(sql.str());
	}
	if (status) {
		std::stringstream sql;
		sql << "CREATE VIEW IF NOT EXISTS loop_overall AS SELECT * FROM parameters p INNER JOIN data_loop_overall d ON p.rowid = d.parameters_fk";
		LOG(trace, "SQL: " << sql.str());
		status &= database->executeSQLStatement(sql.str());
	}
	return status;
}

size_t IterationInformation::prepareParametersTable(
		const size_t _innerSize,
		const size_t _outerSize,
		const MatrixFactorizerOptions &_opts
		)
{
	Table& parameter_table = database->addTable("parameters");
	Table::Tuple_t& parameter_tuple = parameter_table.getTuple();
	parameter_tuple.insert( std::make_pair("k", (int)_innerSize), Table::Parameter);
	parameter_tuple.insert( std::make_pair("n", (int)_outerSize), Table::Parameter);
	parameter_tuple.insert( std::make_pair("p", _opts.px), Table::Parameter);
	parameter_tuple.insert( std::make_pair("r", _opts.py), Table::Parameter);
	parameter_tuple.insert( std::make_pair("sparse_dim", (int)_opts.sparse_dim), Table::Parameter);
	addTupleParametersToTuple(_opts.tuple_parameters, parameter_tuple);
	// we need to add it, otherwise we cannot use checks for presence
	// as they rely on complete table info
	parameter_table.addTuple(parameter_tuple);

	size_t rowid = 0;
	if (database->isDatabaseFileGiven()) {
		// check for presence
		if (!database->isTuplePresentInTable(parameter_table, parameter_tuple)) {
			LOG(debug, "Parameter tuple not present, adding to table.");
			database->writeTable(parameter_table);
		}
		// and return
		rowid = database->getIdOfTuplePresentInTable(
					parameter_table, parameter_tuple);
		// clear table such that present tuple is not stored again
		database->clearTable(parameter_table.getName());
		LOG(info, "Obtaining parameter_key " << rowid);
	} else {
		// else set rowid to arbitrary value as there is no file anyway
		rowid = 1;
	}
	return rowid;
}
