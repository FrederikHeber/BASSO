/*
 * GeneralMinimizer.cpp
 *
 *  Created on: Apr 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "GeneralMinimizer.hpp"

#include "MatrixIO/MatrixIO.hpp"

#include <boost/assign.hpp>
#include <boost/chrono.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/log/trivial.hpp>
#include <cassert>
#include <fenv.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/InverseProblems/QuickAccessReferences.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Norms/Specifics/RegularizedL1Norm.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMappingFactory.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

//#define BREGMANDISTANCEERRORTHRESHOLD 1

using namespace boost::assign;

// static entities
GeneralMinimizer::MinLib_names_t GeneralMinimizer::MinLib_names;

GeneralMinimizer::GeneralMinimizer(
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		const unsigned int _maxinneriter,
		Database &_database,
		const unsigned int _outputsteps
		) :
	Delta(_Delta),
	MaxWalltime(0.),
	MaxOuterIterations(_maxiter),
	MaxInnerIterations(_maxinneriter),
	TolX(1e-6),
	TolY(Delta),
	TolFun(1e-12),
	outputsteps(_outputsteps),
	MinLib(gnuscientificlibrary),
	OldBregmanDistance(0.),
	l2norm(NormFactory::createLpInstance(
			_inverseproblem->x->getSpace(), 2.)),
	parameter_key(0),
	database(_database),
	parameters_table(database.addTable("parameters")),
	data_per_iteration_table(database.addTable("data_per_iteration")),
	data_overall_table(database.addTable("data_overall"))
{
	// set tolerances values
	_inverseproblem->x->getSpace()->getDualityMapping()->setTolerance(TolX);
	_inverseproblem->x->getSpace()->getDualSpace()->getDualityMapping()->setTolerance(TolX);
	_inverseproblem->y->getSpace()->getDualityMapping()->setTolerance(TolY);

	// initalize list with static names
	if (MinLib_names.empty())
		MinLib_names +=
			std::make_pair( "gsl", gnuscientificlibrary),
			std::make_pair( "nlopt", nonlinearoptimization);
}

GeneralMinimizer::~GeneralMinimizer()
{
	// add view to database if not present
	if (!createViews())
		BOOST_LOG_TRIVIAL(warning)
			<< "Could not create overall or per_iteration views.";
}

void GeneralMinimizer::SearchDirection::update(
	 	 const QuickAccessReferences &_refs,
	 	 const SpaceElement_ptr_t &_residual)
{
		Jw = _refs.j_r( _residual );
		BOOST_LOG_TRIVIAL(trace)
			<< "Jw= j_r (R_n) is " << Jw;
		u = _refs.A_t * Jw;
		if (u->getSpace()->getDimension() > 10)
			BOOST_LOG_TRIVIAL(trace)
					<< "newdir is " << u;
		else
			BOOST_LOG_TRIVIAL(debug)
					<< "newdir is " << u;
}

bool GeneralMinimizer::CheckWalltime(
		const boost::chrono::duration<double> &_time) const
{
	if (MaxWalltime != boost::chrono::duration<double>(0.))
		return _time >= MaxWalltime;
	else
		return false;
}

bool GeneralMinimizer::CheckIterations(
		const int _current_outeriterations) const
{
	// walltime overrules maxiter
	if (MaxWalltime == boost::chrono::duration<double>(0.))
		return _current_outeriterations >= MaxOuterIterations;
	else
		return false;
}

bool GeneralMinimizer::CheckResiduum(
		const double _residuum) const
{
	return _residuum <= TolY;
}

bool GeneralMinimizer::CheckRelativeResiduum(
		const double _residuum,
		const double _ynorm) const
{
	return _residuum/_ynorm <= TolY;
}

void GeneralMinimizer::ReturnValues::output(
		const double ynorm) const
{
	/// output prior to iterate update
	BOOST_LOG_TRIVIAL(debug)<< "#" << NumberOuterIterations
	<< " with residual of " << residuum;
	BOOST_LOG_TRIVIAL(debug)
	<< "#" << NumberOuterIterations << ": "
	<< "||Ax_n-y||/||y|| is " << residuum/ynorm;
	BOOST_LOG_TRIVIAL(trace)
	<< "x_n is " << m_solution;
	BOOST_LOG_TRIVIAL(trace)
	<< "R_n is " << m_residual;
}

double GeneralMinimizer::calculateResidual(
		const InverseProblem_ptr_t &_problem,
		SpaceElement_ptr_t &_residual
		) const
{
	const LinearMapping &A = static_cast<const LinearMapping &>(*_problem->A);
	*_residual = A * _problem->x;
	*_residual -= _problem->y;
	const Norm &NormY = *_problem->y->getSpace()->getNorm();
	return NormY(_residual);
}

const double GeneralMinimizer::calculateBregmanDistance(
		const boost::shared_ptr<BregmanDistance> &_Delta_p,
		const SpaceElement_ptr_t &_solution,
		const SpaceElement_ptr_t &_truesolution,
		const SpaceElement_ptr_t &_dual_solution) const
{
	double distance = 0.;
	if (!_truesolution->isZero()) {
		distance = (*_Delta_p)(
			_solution,
			_truesolution,
			_dual_solution);
#ifdef BREGMANDISTANCEERRORTHRESHOLD
		const int roundmode = fegetround();
		fesetround(FE_DOWNWARD);
		const double lower_bound = (*_Delta_p)(
			_solution,
			_truesolution,
			_dual_solution);
		fesetround(FE_UPWARD);
		const double upper_bound = (*_Delta_p)(
			_solution,
			_truesolution,
			_dual_solution);
		fesetround(roundmode);
		const double errorvalue =
				std::max(distance-lower_bound, upper_bound-distance);
		BOOST_LOG_TRIVIAL(debug)
				<< "Reduction in Bregman Distance is " << OldBregmanDistance-distance;
		BOOST_LOG_TRIVIAL(debug)
				<< "Bregman distance is " << distance
				<< "+-" << errorvalue;
//				<< " in [" << lower_bound << "," << upper_bound << "]";
#else
		BOOST_LOG_TRIVIAL(debug)
				<< "Bregman distance is " << distance;
#endif
		// check that distance truly decreases
		assert( (OldBregmanDistance == 0.)
				|| ((OldBregmanDistance - distance) > - BASSOTOLERANCE) );
		OldBregmanDistance = distance;
	}
	return distance;
}
const double GeneralMinimizer::calculateError(
		const SpaceElement_ptr_t &_solution,
		const SpaceElement_ptr_t &_truesolution) const
{
	const Norm &NormX = *_solution->getSpace()->getNorm();
	double new_error = 0.;
	if (!_truesolution->isZero()) {
		if (static_cast<const RegularizedL1Norm *>(&NormX) == NULL) {
			new_error = NormX(_solution-_truesolution);
		} else {
			// create L2 norm for measuring error
			new_error = (*l2norm)(_solution-_truesolution);
		}
		BOOST_LOG_TRIVIAL(debug)
			<< "Error is " << ": "
			<< "||x_n-x|| is " << new_error;
	}
	return new_error;
}


bool GeneralMinimizer::isValidMinLibName(const std::string &_name)
{
	MinLib_names_t::const_iterator iter =
			MinLib_names.find(_name);
	return (iter != MinLib_names.end());
}

void GeneralMinimizer::setMinLib(const std::string &_name)
{
	assert( isValidMinLibName(_name) );
	MinLib = MinLib_names[_name];
}

void GeneralMinimizer::printIntermediateSolution(
		const SpaceElement_ptr_t &_solution,
		const LinearMapping &_A,
		unsigned int _NumberOuterIterations
		) const
{
	// print each solution
	if ((outputsteps != 0) &&
			(_NumberOuterIterations % outputsteps == 0)) {
		{
			std::stringstream solution_file;
			solution_file << "solution"
					<< (_NumberOuterIterations / outputsteps) << ".m";
			using namespace MatrixIO;
			std::ofstream ost(solution_file.str().c_str());
			if (ost.good())
				try {
					ost << _solution;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Could not write all data of intermediate solution to stream.\n";
				}
			else {
				std::cerr << "Failed to open " << solution_file.str() << std::endl;
			}
		}
		{
			std::stringstream solution_file;
			solution_file << "projected_solution"
					<< (_NumberOuterIterations / outputsteps) << ".m";
			using namespace MatrixIO;
			std::ofstream ost(solution_file.str().c_str());
			if (ost.good())
				try {
					ost << _A * _solution;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Could not write all data of projected intermediate solution to stream.\n";
				}
			else {
				std::cerr << "Failed to open " << solution_file.str() << std::endl;
			}
		}

	}
}

void GeneralMinimizer::setParameterKey(
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
		addAdditionalParametersToTuple(parameter_tuple, do_replace);
		// add additional parameters specified from user
		for (std::vector<std::string>::const_iterator iter = tuple_params.begin();
				iter != tuple_params.end(); ) {
			const std::string &token = (*iter++);
			const std::string &value = (*iter++);
			try {
				const int int_value = boost::lexical_cast<int>(value);
				const double double_value = boost::lexical_cast<double>(value);
				if ((double)int_value == double_value)
					parameter_tuple.insert( std::make_pair(token, int_value), Table::Parameter);
				else
					parameter_tuple.insert( std::make_pair(token, double_value), Table::Parameter);
			} catch(const boost::bad_lexical_cast &) {
				parameter_tuple.insert( std::make_pair(token, value), Table::Parameter);
			}
		}
	} else {
		parameter_tuple.replace( "p", _val_NormX);
		parameter_tuple.replace( "r", _val_NormY);
		parameter_tuple.replace( "N", (int)_N);
		parameter_tuple.replace( "dim", (int)_dim);
		parameter_tuple.replace( "max_iterations", _MaxOuterIterations);
		// add additional parameters from derived classes
		addAdditionalParametersToTuple(parameter_tuple, do_replace);
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
			BOOST_LOG_TRIVIAL(debug)
					<< "Parameter tuple not present, adding to table.";
			database.writeTable(parameters_table);
		}
		// and store rowid
		rowid = database.getIdOfTuplePresentInTable(
					parameters_table, parameter_tuple);
		// clear table such that present tuple is not stored again
		database.clearTable(parameters_table.getName());
		BOOST_LOG_TRIVIAL(info)
			<< "Setting parameter_key to " << rowid;
	} else {
		// else set rowid to arbitrary value as there is no file anyway
		rowid = 1;
	}
	const_cast<size_t &>(parameter_key) = rowid;
}

Table::Tuple_t & GeneralMinimizer::preparePerIterationTuple() const
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

Table::Tuple_t & GeneralMinimizer::prepareOverallTuple() const
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

void GeneralMinimizer::finalizeOverallTuple(
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

bool GeneralMinimizer::createViews() const
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
	if (!status)
		BOOST_LOG_TRIVIAL(error)
			<< "(Some of the) Required Tables are empty, not creating views.";
	if (status) {
		std::stringstream sql;
		sql << "CREATE VIEW IF NOT EXISTS overall AS SELECT * FROM parameters p INNER JOIN data_overall d ON p.rowid = d.parameters_fk";
		BOOST_LOG_TRIVIAL(trace)
			<< "SQL: " << sql.str();
		status &= database.executeSQLStatement(sql.str());
	}
	if (status) {
		std::stringstream sql;
		sql << "CREATE VIEW IF NOT EXISTS per_iteration AS SELECT * FROM parameters p INNER JOIN data_per_iteration d ON p.rowid = d.parameters_fk";
		BOOST_LOG_TRIVIAL(trace)
			<< "SQL: " << sql.str();
		status &= database.executeSQLStatement(sql.str());
	}
	return status;
}

void GeneralMinimizer::setAdditionalTupleParameters(
			const std::vector<std::string> &_tuple_params)
{
	assert ( _tuple_params.size() % 2 == 0);
	const_cast< std::vector<std::string> &>(tuple_params) = _tuple_params;
}

void GeneralMinimizer::resetState()
{
	// reset Bregman distance check variable
	resetBregmanDistance();
	// reset state of derived objects
	resetState_interal();
}
