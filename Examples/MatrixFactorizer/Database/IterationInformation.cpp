/*
 * IterationInformation.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "IterationInformation.hpp"

#include <boost/lexical_cast.hpp>

#include "Log/Logging.hpp"
#include "MatrixFactorizer/Options/MatrixFactorizerOptions.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"

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
	// install parameter set and get its key
	const size_t parameter_key = prepareParametersTable(
			_innersize,
			_outersize,
			_opts);

	loop_tuple.insert( std::make_pair("parameters_fk", (int)parameter_key), Table::Parameter);
	loop_tuple.insert( std::make_pair("loop_nr", (int)0), Table::Data);
	loop_tuple.insert( std::make_pair("residual", 0.), Table::Data);
	loop_tuple.insert( std::make_pair("scaling_change", 0.), Table::Data);

	overall_tuple.insert( std::make_pair("parameters_fk", (int)parameter_key), Table::Parameter);
	overall_tuple.insert( std::make_pair("loops", (int)0), Table::Data);
	overall_tuple.insert( std::make_pair("residual", 0.), Table::Data);
	overall_tuple.insert( std::make_pair("runtime", 0.), Table::Data);
}

IterationInformation::~IterationInformation()
{
	// add views after tables have been created and filled
	if (!createViews())
		BOOST_LOG_TRIVIAL(warning)
			<< "Could not create overall or per_iteration views.";
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
		BOOST_LOG_TRIVIAL(error)
			<< "Unknown table in IterationInformation::addTuple()";
	}
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
	if (!status)
		BOOST_LOG_TRIVIAL(error)
			<< "(Some of the) Required Tables are empty, not creating views.";
	if (status) {
		std::stringstream sql;
		sql << "CREATE VIEW IF NOT EXISTS loop AS SELECT * FROM parameters p INNER JOIN data_loop d ON p.rowid = d.parameters_fk";
		BOOST_LOG_TRIVIAL(trace)
			<< "SQL: " << sql.str();
		status &= database->executeSQLStatement(sql.str());
	}
	if (status) {
		std::stringstream sql;
		sql << "CREATE VIEW IF NOT EXISTS loop_overall AS SELECT * FROM parameters p INNER JOIN data_loop_overall d ON p.rowid = d.parameters_fk";
		BOOST_LOG_TRIVIAL(trace)
			<< "SQL: " << sql.str();
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
	for (std::vector<std::string>::const_iterator iter = _opts.tuple_parameters.begin();
			iter != _opts.tuple_parameters.end(); ) {
		const std::string &token = (*iter++);
		const std::string &value = (*iter++);
		try {
			const int int_value = boost::lexical_cast<int>(value);
			parameter_tuple.insert( std::make_pair(token, int_value), Table::Parameter);
			BOOST_LOG_TRIVIAL(debug)
					<< " Adding additional integer parameter ("
					<< token << "," << int_value << ") to loop tuple.";
		} catch(const boost::bad_lexical_cast &) {
			try {
				const double double_value = boost::lexical_cast<double>(value);
				parameter_tuple.insert( std::make_pair(token, double_value), Table::Parameter);
				BOOST_LOG_TRIVIAL(debug)
						<< " Adding additional double parameter ("
						<< token << "," << double_value << ") to loop tuple.";
			} catch(const boost::bad_lexical_cast &) {
				parameter_tuple.insert( std::make_pair(token, value), Table::Parameter);
				BOOST_LOG_TRIVIAL(debug)
						<< " Adding additional string parameter ("
						<< token << "," << value << ") to loop tuple.";
			}
		}
	}
	// we need to add it, otherwise we cannot use checks for presence
	// as they rely on complete table info
	parameter_table.addTuple(parameter_tuple);

	size_t rowid = 0;
	if (database->isDatabaseFileGiven()) {
		// check for presence
		if (!database->isTuplePresentInTable(parameter_table, parameter_tuple)) {
			BOOST_LOG_TRIVIAL(debug)
					<< "Parameter tuple not present, adding to table.";
			database->writeTable(parameter_table);
		}
		// and return
		rowid = database->getIdOfTuplePresentInTable(
					parameter_table, parameter_tuple);
		// clear table such that present tuple is not stored again
		database->clearTable(parameter_table.getName());
		BOOST_LOG_TRIVIAL(info)
			<< "Obtaining parameter_key " << rowid;
	} else {
		// else set rowid to arbitrary value as there is no file anyway
		rowid = 1;
	}
	return rowid;
}
