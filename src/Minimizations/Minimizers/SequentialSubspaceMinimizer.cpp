/*
 * SequentialSubspaceMinimizer.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "SequentialSubspaceMinimizer.hpp"

#include <boost/chrono.hpp>
#include <boost/assign.hpp>
#include <boost/log/trivial.hpp>
#include <cassert>
#include <cmath>
#include <iterator>
#include <numeric>
#include <sstream>

#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/InverseProblems/QuickAccessReferences.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/HyperplaneProjection.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizerFactory.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Functions/Minimizers/MinimizerExceptions.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Options/CommandLineOptions.hpp"

using namespace boost::assign;

SequentialSubspaceMinimizer::SequentialSubspaceMinimizer(
		const CommandLineOptions &_opts,
		const InverseProblem_ptr_t &_inverseproblem,
		Database &_database
		) :
	GeneralMinimizer(
			_opts,
			_inverseproblem,
			_database
			),
	N(2),
	inexactLinesearch(false),
	constant_positivity(1e-6),
	constant_interpolation(0.6),
	DoCalculateAngles(false),
	OrthogonalizationType(_opts.orthogonalization_type)
{}

void SequentialSubspaceMinimizer::setN(
		const unsigned int _N
		)
{
	// check that number of search directions is greater than 1
	if (_N < 1)
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("N");
	const_cast<unsigned int&>(N) = _N;
	// make state invalid
	istate.reset();
}

static Table::Tuple_t& prepareAngleTuple(
		Table &_table,
		const int _parameter_key,
		const size_t _maxangles)
{
	assert(_parameter_key != 0);

	Table::Tuple_t &angle_tuple = _table.getTuple();
	angle_tuple.insert( std::make_pair("parameters_fk", (int)_parameter_key), Table::Parameter);
	angle_tuple.insert( std::make_pair("iteration", (int)0), Table::Data);
	std::vector<std::string> names;
	names += "angle","bregman_angle";
	for (std::vector<std::string>::const_iterator nameiter = names.begin();
			nameiter != names.end(); ++nameiter)
		for (unsigned int i=0; i<_maxangles; ++i) {
			std::stringstream componentname;
			componentname << *nameiter << i+1;
			BOOST_LOG_TRIVIAL(debug)
				<< "Adding " << componentname.str();
			angle_tuple.insert( std::make_pair(componentname.str(), 0.), Table::Data);
		}
	return angle_tuple;
}


bool createAnglesViews(const Database &_database)
{
	// write tables beforehand
	_database.writeAllTables();
	// we create views angles (which were present before the switch to distinct
	// data and parameters table)
	bool status = true;
	{
		// check whether tables are present and contain elements
		const Table &data_angles_table = _database.getTableConst("data_angles");
		status &= !data_angles_table.empty();
		// we don't  check for parameter table's non-emptiness as is
		// possibly might be if the used parameter tuple is already in
		// the database, see setParameterKey()
	}
	if (!status)
		BOOST_LOG_TRIVIAL(error)
			<< "(Some of the) Required Tables are empty, not creating angle views.";
	if (status) {
		std::stringstream sql;
		sql << "CREATE VIEW IF NOT EXISTS angles AS SELECT * FROM parameters p INNER JOIN data_angles d ON p.rowid = d.parameters_fk";
		BOOST_LOG_TRIVIAL(trace)
			<< "SQL: " << sql.str();
		status &= _database.executeSQLStatement(sql.str());
	}
	return status;
}

void
SequentialSubspaceMinimizer::addAdditionalParametersToTuple(
		Table::Tuple_t &_tuple,
		const bool _do_replace) const
{
	if (!_do_replace) {
		if (inexactLinesearch) {
			_tuple.insert( std::make_pair("c1", constant_positivity), Table::Parameter);
			_tuple.insert( std::make_pair("c2", constant_interpolation), Table::Parameter);
		}
		_tuple.insert( std::make_pair("max_inner_iterations", MaxInnerIterations), Table::Parameter);
	} else {
		if (inexactLinesearch) {
			_tuple.replace( "c1", constant_positivity);
			_tuple.replace( "c2", constant_interpolation);
		}
		_tuple.replace( "max_inner_iterations", MaxInnerIterations);
	}
}

bool SequentialSubspaceMinimizer::isNonConverging(
		const double current_residuum,
		const double initial_residuum) const
{
	static const double threshold_factor = 1e4;
	/// check for non-convergence
	if (current_residuum/initial_residuum >= threshold_factor) {
		BOOST_LOG_TRIVIAL(info)<< "Current residuum is "
		<< current_residuum
		<< " exceeding initial value of 1. by "
		<< current_residuum/initial_residuum;
		BOOST_LOG_TRIVIAL(error)
		<< "STOPPING ITERATION";
		return true;
	} else
		return false;
}

void SequentialSubspaceMinimizer::fillPerIterationTable(
		Table::Tuple_t& per_iteration_tuple)
{
	// create a sequence of entries in Database till MaxOuterIterations
	for (;istate.NumberOuterIterations < MaxOuterIterations;
			++istate.NumberOuterIterations) {
		// update iterations in tuple
		per_iteration_tuple.replace( "iteration", (int)istate.NumberOuterIterations);

		// submit current tuples
		data_per_iteration_table.addTuple(per_iteration_tuple);
	}
}

void SequentialSubspaceMinimizer::fillAngleTable(
		Table::Tuple_t& angle_tuple,
		Table& data_angle_table)
{
	// create a sequence of entries in Database till MaxOuterIterations
	for (;istate.NumberOuterIterations < MaxOuterIterations;
			++istate.NumberOuterIterations) {
		// update iterations in tuple
		if (DoCalculateAngles)
		angle_tuple.replace( "iteration", (int)istate.NumberOuterIterations);

		// submit current tuples
		if (DoCalculateAngles)
			data_angle_table.addTuple(angle_tuple);
	}
}

const unsigned int SequentialSubspaceMinimizer::calculateStepWidth(
		const QuickAccessReferences& refs,
		const SpaceElement_ptr_t& dual_solution,
		std::vector<double> & tmin,
		const std::vector<SpaceElement_ptr_t> &_searchspace,
		const std::vector<double> &_alphas
		) const
{
	// consistency checks
	const size_t Ndirs = _searchspace.size();
	assert( Ndirs == _alphas.size() );
	assert( Ndirs == tmin.size() );

	// tmin=fminunc(@(t) BregmanProjectionFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
	BregmanProjectionFunctional bregman(
			refs.DualNormX, refs.J_q, refs.J_q.getPower(), _searchspace, _alphas);
	const HyperplaneProjection<BregmanProjectionFunctional> functional(
			bregman, dual_solution);

	FunctionalMinimizer< std::vector<double> >::ptr_t functionminimizer;
	if (inexactLinesearch) {
		const unsigned int current_index = istate.searchspace->getIndex();
		Wolfe_indexset_t Wolfe_indexset(1, current_index);
		functionminimizer =
				FunctionalMinimizerFactory::create< std::vector<double> >(
					Ndirs,
					functional,
					tmin,
					constant_positivity,
					Wolfe_indexset,
					constant_interpolation);
	} else {
		functionminimizer =
				FunctionalMinimizerFactory::create< std::vector<double> >(
					Ndirs,
					functional,
					tmin);
	}
	functionminimizer->setMaxIterations(MaxInnerIterations);
	unsigned int inner_iterations = (*functionminimizer)(
		Ndirs, TolFun, tmin);

	std::stringstream output_stepwidth;
	std::copy(tmin.begin(), tmin.end(),
			std::ostream_iterator<double>(output_stepwidth, " "));
	BOOST_LOG_TRIVIAL(debug)<< "tmin " << output_stepwidth.str()
	<< " found in " << inner_iterations << " iterations.";

	return inner_iterations;
}

void SequentialSubspaceMinimizer::updateAngleTable(
		const SpaceElement_ptr_t& newdir,
		Table::Tuple_t& angle_tuple) const
{
	if (DoCalculateAngles) {
		const LastNSearchDirections& searchspace =
				static_cast<const LastNSearchDirections &>(*istate.searchspace);
		const std::vector<unsigned int> &lastIndices = searchspace.getLastIndices();
		angle_tuple.replace("iteration", (int) (istate.NumberOuterIterations));
		// calculate bregman angles for angles database
		{
			const IterationState::angles_t angles =
					istate.calculateBregmanAngles(newdir);
			for (unsigned int i = 0; i < angles.size();
					++i) {
				std::stringstream componentname;
				componentname << "bregman_angle" << lastIndices[i] + 1;
				angle_tuple.replace(componentname.str(), angles[i]);
			}
		}
		// calculate "scalar product" angles for angles database
		{
			const IterationState::angles_t angles = istate.calculateAngles(
					newdir);
			for (unsigned int i = 0; i < angles.size();
					++i) {
				std::stringstream componentname;
				componentname << "angle" << lastIndices[i] + 1;
				angle_tuple.replace(componentname.str(), angles[i]);
			}
		}
	}
}

void SequentialSubspaceMinimizer::updateSearchspace(
		const SpaceElement_ptr_t& _truesolution,
		const SpaceElement_ptr_t& newdir,
		const double alpha)
{
	/// update search space with new direction
	if (_truesolution->isZero()) {
		istate.updateSearchSpace(newdir, alpha, istate.m_dual_solution,
				istate.m_solution);
	} else {
		istate.updateSearchSpace(newdir, alpha, istate.m_dual_solution, _truesolution);
	}
	BOOST_LOG_TRIVIAL(trace)<< "updated_index is " << istate.searchspace->getIndex();
}

void SequentialSubspaceMinimizer::updateIterates(
		const QuickAccessReferences& refs,
		const std::vector<double> tmin,
		SpaceElement_ptr_t& _x,
		SpaceElement_ptr_t& dual_x) const
{
	// x=DualityMapping(Jx-tmin*u,DualNormX,DualPowerX,TolX);
	{
		const SpaceElement_ptr_t tempelement = refs.DualSpaceX.createElement();
		for (size_t i = 0; i < N; ++i)
			*tempelement += tmin[i] * istate.getSearchSpace()[i];
		*istate.m_dual_solution -= tempelement;
	}
	BOOST_LOG_TRIVIAL(trace)<< "x^*_n+1 is " << istate.m_dual_solution;
	*istate.m_solution = refs.J_q(istate.m_dual_solution);
	BOOST_LOG_TRIVIAL(trace)<< "x_n+1 is " << istate.m_solution;
	*_x = istate.m_solution;
	*dual_x = istate.m_dual_solution;
}

SequentialSubspaceMinimizer::ReturnValues
SequentialSubspaceMinimizer::operator()(
		const InverseProblem_ptr_t &_problem,
		const SpaceElement_ptr_t &_startvalue,
		const SpaceElement_ptr_t &_dualstartvalue,
		const SpaceElement_ptr_t &_truesolution
		)
{
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	// gather some refs for convenient access
	QuickAccessReferences refs(_problem);

	/// -# initialize return structure
	if (!istate.getisInitialized()) {
		SpaceElement_ptr_t residual = refs.SpaceY.createElement();
		const double residuum = calculateResidual(
			_problem,
			residual);
		istate.set(
				_startvalue,
				_dualstartvalue,
				residual,
				residuum,
				N,
				OrthogonalizationType);
	}

	// create Bregman distance object
	boost::shared_ptr<BregmanDistance> Delta_p;
	if (!_truesolution->isZero())
		Delta_p.reset(new BregmanDistance (
				refs.NormX, refs.J_p, refs.J_p.getPower()));

	/// build data tuple for iteration, overall, and angles information
	setParameterKey(
			refs.NormX.getPvalue(),
			refs.NormY.getPvalue(),
			N,
			refs.SpaceX.getDimension(),
			MaxOuterIterations);
	Table::Tuple_t& per_iteration_tuple = preparePerIterationTuple(
			data_per_iteration_table, parameter_key);
	per_iteration_tuple.insert( std::make_pair("inner_iterations", (int)0), Table::Data);
	Table::Tuple_t& overall_tuple = prepareOverallTuple(
			data_overall_table, parameter_key);
	Table& data_angle_table = database.addTable("data_angles");
	Table::Tuple_t& angle_tuple = prepareAngleTuple(
			data_angle_table,
			parameter_key,
			N);

	/// -# check initial stopping criterion
	const double ynorm = refs.NormY(refs.y);
	const double initial_relative_residuum = fabs(istate.residuum/ynorm);
	bool StopCriterion = CheckStoppingCondition(
			boost::chrono::duration<double>(0.),
			istate.NumberOuterIterations,
			istate.residuum,
			ynorm);

	istate.status = ReturnValues::started;
	while (!StopCriterion) {
		/// Calculation of search direction
		// Jw=DualityMapping(w,NormY,PowerY,TolX);
		searchdir.update(refs, istate.m_residual);
		if (searchdir.u->isZero(BASSOTOLERANCE)) {
			StopCriterion = true;
			BOOST_LOG_TRIVIAL(debug)
				<< "newdir is zero, stopping.";
			continue;
		}
		updateAngleTable(searchdir.u, angle_tuple);

		/// output prior to iterate update
		istate.output(ynorm);

		/// update search space with new direction
		// alpha=Jw'*y
		const double alpha =
				searchdir.Jw * refs.y;
		BOOST_LOG_TRIVIAL(trace)
			<< "alpha is " << alpha;

		// add u to U and alpha to alphas
		updateSearchspace(_truesolution, searchdir.u, alpha);

		/// get optimal stepwidth
		std::vector<double> tmin(N, 0.);
		double stepwidth_norm = 0.;
		unsigned int inner_iterations = 0;
		try {
			inner_iterations =
					calculateStepWidth(refs, istate.m_dual_solution, tmin,
							istate.getSearchSpace(), istate.getAlphas());
			stepwidth_norm = std::inner_product(tmin.begin(), tmin.end(), tmin.begin(), stepwidth_norm);
		} catch (MinimizerIllegalNumber_exception &e) {
			BOOST_LOG_TRIVIAL(error)
					<< "Encountered illegal number in line search minimum, not updating.";
			tmin = std::vector<double>(N,0.);
		}

#ifdef FULLMATRIXNORM
		/// give estimate on reduction in Bregman distance
		if (!_truesolution->isZero()) {
			const double C = 1.;
			const double p = refs.NormX.getPvalue();
			const double reduction =
					::pow(istate.residuum, p)
					 /(p*::pow(C,p-1.)*::pow(refs.A.Norm(),p));
			BOOST_LOG_TRIVIAL(debug)
					<< "Estimated reduction in Bregman distance is " << reduction;
		}
#endif

		/// database update prior to iterate update
		per_iteration_tuple.replace( "iteration", (int)istate.NumberOuterIterations);
		per_iteration_tuple.replace( "residual", istate.residuum);
		per_iteration_tuple.replace( "relative_residual", istate.residuum/ynorm);
		per_iteration_tuple.replace( "bregman_distance",
				calculateBregmanDistance(
								Delta_p, istate.m_solution, _truesolution, istate.m_dual_solution));
		per_iteration_tuple.replace( "error",
				calculateError(istate.m_solution, _truesolution));
		per_iteration_tuple.replace( "updated_index", (int)istate.searchspace->getIndex());
		per_iteration_tuple.replace("inner_iterations",
				(int) (inner_iterations));
		per_iteration_tuple.replace( "stepwidth", sqrt(stepwidth_norm));
		// submit current tuples
		data_per_iteration_table.addTuple(per_iteration_tuple);
		if (DoCalculateAngles)
			data_angle_table.addTuple(angle_tuple);

		/// update iterate
		updateIterates(refs, tmin, _problem->x, istate.m_dual_solution);

		/// update residual
		istate.residuum = calculateResidual(_problem,istate.m_residual);

		/// check iterations count/wall time
		boost::chrono::high_resolution_clock::time_point timing_intermediate =
				boost::chrono::high_resolution_clock::now();
		++istate.NumberOuterIterations;
		StopCriterion = CheckStoppingCondition(
			timing_intermediate - timing_start,
			istate.NumberOuterIterations,
			istate.residuum,
			ynorm);

		/// check for non-convergence
		const double current_relative_residuum = fabs(istate.residuum/ynorm);
		if (isNonConverging(current_relative_residuum,
				initial_relative_residuum)) {
			fillPerIterationTable(per_iteration_tuple);
			fillAngleTable(angle_tuple,data_angle_table);
		}

		// print intermediate solution
		printIntermediateSolution(
				istate.m_solution, refs.A, istate.NumberOuterIterations);
	}
	/// last output of istate
	istate.output(ynorm);

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();

	// submit overall_tuple
	overall_tuple.replace( "iterations", istate.NumberOuterIterations );
	overall_tuple.replace( "residual", istate.residuum );
	overall_tuple.replace( "relative_residual", istate.residuum/ynorm );
	overall_tuple.replace( "runtime",
			boost::chrono::duration<double>(timing_end - timing_start).count() );
	finalizeOverallTuple(overall_tuple, refs);
	data_overall_table.addTuple(overall_tuple);

	// create angles view if desired
	if ((DoCalculateAngles) && (!createAnglesViews(database)))
		BOOST_LOG_TRIVIAL(warning)
			<< "Could not create angles view in SQLite database.";

	// and return solution
	istate.status = ReturnValues::finished;
	return static_cast<ReturnValues &>(istate);
}

void
SequentialSubspaceMinimizer::setupdateIndexAlgorithm(
		enum LastNSearchDirections::UpdateAlgorithmType _type)
{
	LastNSearchDirections::updateIndexType = _type;
}

void
SequentialSubspaceMinimizer::setEnforceRandomMapping(
		const bool _enforceRandomMapping)
{
	LastNSearchDirections::enforceRandomMapping = _enforceRandomMapping;
}

