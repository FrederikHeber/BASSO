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
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Functions/Minimizers/Minimizer.hpp"
#include "Minimizations/Functions/Minimizers/MinimizerExceptions.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Norms/RegularizedL1Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

using namespace boost::assign;

SequentialSubspaceMinimizer::SequentialSubspaceMinimizer(
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		const unsigned int _maxinneriter,
		Database &_database,
		const unsigned int _outputsteps
		) :
	GeneralMinimizer(
			_inverseproblem,
			_Delta,
			_maxiter,
			_maxinneriter,
			_database,
			_outputsteps
			),
	N(2),
	inexactLinesearch(false),
	constant_positivity(1e-6),
	constant_interpolation(0.6),
	DoCalculateAngles(false)
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

//!> limit for angle columns in angles table
static const unsigned int MAXANGLES = 4;

static Table::Tuple_t prepareAngleTuple(
		const double _val_NormX,
		const double _val_NormY,
		const unsigned int _N,
		const unsigned int _dim)
{
	Table::Tuple_t angle_tuple;
	angle_tuple.insert( std::make_pair("p", _val_NormX), Table::Parameter);
	angle_tuple.insert( std::make_pair("r", _val_NormY), Table::Parameter);
	angle_tuple.insert( std::make_pair("N", (int)_N), Table::Parameter);
	angle_tuple.insert( std::make_pair("dim", (int)_dim), Table::Parameter);
	for (std::vector<std::string>::const_iterator iter
			= GeneralMinimizer::tuple_params.begin();
			iter != GeneralMinimizer::tuple_params.end(); ++iter) {
		const std::string &token = *iter++;
		const std::string &value = *iter++;
		angle_tuple.insert( std::make_pair(token, value), Table::Parameter);
	}
	angle_tuple.insert( std::make_pair("iteration", (int)0), Table::Data);
	std::vector<std::string> names;
	names += "angle","bregman_angle";
	for (std::vector<std::string>::const_iterator nameiter = names.begin();
			nameiter != names.end(); ++nameiter)
		for (unsigned int i=0; i<MAXANGLES; ++i) {
			std::stringstream componentname;
			componentname << *nameiter << i+1;
			BOOST_LOG_TRIVIAL(debug)
				<< "Adding " << componentname.str();
			angle_tuple.insert( std::make_pair(componentname.str(), 0.), Table::Data);
		}
	return angle_tuple;
}

SpaceElement_ptr_t calculateDualStartingValue(
		const SpaceElement_ptr_t &_dualstartvalue)
{
	SpaceElement_ptr_t dual_solution =
			_dualstartvalue->getSpace()->createElement();
	*dual_solution = _dualstartvalue;
	BOOST_LOG_TRIVIAL(trace)
		<< "Jx_0 is " << dual_solution;
	return dual_solution;
}

Table::Tuple_t SequentialSubspaceMinimizer::addInfoToPerIterationTable(
		const QuickAccessReferences &_refs) const
{
	Table::Tuple_t per_iteration_tuple = preparePerIterationTuple(
			_refs.NormX.getPvalue(),
			_refs.NormY.getPvalue(),
			N,
			_refs.SpaceX.getDimension(),
			MaxOuterIterations);
	if (inexactLinesearch) {
		per_iteration_tuple.insert( std::make_pair("c1", constant_positivity), Table::Parameter);
		per_iteration_tuple.insert( std::make_pair("c2", constant_interpolation), Table::Parameter);
	}
	per_iteration_tuple.insert( std::make_pair("max_inner_iterations", MaxInnerIterations), Table::Parameter);
	per_iteration_tuple.insert( std::make_pair("inner_iterations", (int)0), Table::Data);

	return per_iteration_tuple;
}

Table::Tuple_t SequentialSubspaceMinimizer::addInfoToOverallTable(
		const QuickAccessReferences &_refs) const
{
	Table::Tuple_t overall_tuple = prepareOverallTuple(
			_refs.NormX.getPvalue(),
			_refs.NormY.getPvalue(),
			N,
			_refs.SpaceX.getDimension(),
			MaxOuterIterations);
	if (inexactLinesearch) {
		overall_tuple.insert( std::make_pair("c1", constant_positivity), Table::Parameter);
		overall_tuple.insert( std::make_pair("c2", constant_interpolation), Table::Parameter);
	}
	overall_tuple.insert( std::make_pair("max_inner_iterations", MaxInnerIterations), Table::Parameter);
	// due to Eigen's lazy evaluation runtime is not measured accurately

	return overall_tuple;
}

Table::Tuple_t SequentialSubspaceMinimizer::addInfoToAnglesTable(
		const QuickAccessReferences &_refs) const
{
	Table::Tuple_t angle_tuple;
	if (DoCalculateAngles) {
		// build angle tuple for search direction angle information
		angle_tuple = prepareAngleTuple(
				_refs.NormX.getPvalue(),
				_refs.NormY.getPvalue(),
				N,
				_refs.SpaceX.getDimension());
		angle_tuple.insert( std::make_pair("max_iterations", MaxOuterIterations), Table::Parameter);
		angle_tuple.insert( std::make_pair("max_inner_iterations", MaxInnerIterations), Table::Parameter);
	}
	return angle_tuple;
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
		Table::Tuple_t& per_iteration_tuple,
		Table& per_iteration_table)
{
	// create a sequence of entries in Database till MaxOuterIterations
	for (;istate.NumberOuterIterations < MaxOuterIterations;
			++istate.NumberOuterIterations) {
		// update iterations in tuple
		per_iteration_tuple.replace( "iteration", (int)istate.NumberOuterIterations);

		// submit current tuples
		per_iteration_table.addTuple(per_iteration_tuple);
	}
}

void SequentialSubspaceMinimizer::fillAngleTable(
		Table::Tuple_t& angle_tuple,
		Table& angle_table)
{
	// create a sequence of entries in Database till MaxOuterIterations
	for (;istate.NumberOuterIterations < MaxOuterIterations;
			++istate.NumberOuterIterations) {
		// update iterations in tuple
		if (DoCalculateAngles)
		angle_tuple.replace( "iteration", (int)istate.NumberOuterIterations);

		// submit current tuples
		if (DoCalculateAngles)
			angle_table.addTuple(angle_tuple);
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
	BregmanProjectionFunctional bregman(refs.DualNormX,
			dynamic_cast<const PowerTypeDualityMapping&>(refs.J_q),
			refs.J_q.getPower());
	const HyperplaneProjection functional(bregman, dual_solution,
			_searchspace, _alphas);

	// due to templation we need to instantiate both, as user
	// decides during runtime which we need
	Minimizer<gsl_vector> minimizer_gsl(Ndirs);
	Minimizer<NLopt_vector> minimizer_nlopt(Ndirs);
	const unsigned int current_index = istate.searchspace->getIndex();
	Wolfe_indexset_t Wolfe_indexset(1, current_index);
	if (MinLib == gnuscientificlibrary) {
		minimizer_gsl.setMaxIterations(MaxInnerIterations);
	} else if (MinLib == nonlinearoptimization) {
		// hand index set to minimizer
		minimizer_nlopt.setConstantPositivity(constant_positivity);
		minimizer_nlopt.setPositivityBoundIndices(Wolfe_indexset);
		minimizer_nlopt.setMaxIterations(MaxInnerIterations);
	}

	unsigned int inner_iterations = 0;
	if (inexactLinesearch) {
		FunctionalMinimizer<std::vector<double>, gsl_vector> fmin_gsl(
				functional, minimizer_gsl, tmin, constant_positivity,
				constant_interpolation);
		FunctionalMinimizer<std::vector<double>, NLopt_vector> fmin_nlopt(
				functional, minimizer_nlopt, tmin, constant_positivity,
				constant_interpolation);
		switch (MinLib) {
		case gnuscientificlibrary:
			inner_iterations = fmin_gsl(Ndirs, TolFun, Wolfe_indexset, tmin);
			break;
		case nonlinearoptimization:
			inner_iterations = fmin_nlopt(Ndirs, TolFun, Wolfe_indexset, tmin);
			break;
		default:
			throw MinimizationIllegalValue_exception()
					<< MinimizationIllegalValue_name("MinLib");
			break;
		}
	} else {
		FunctionalMinimizer<std::vector<double>, gsl_vector> fmin_gsl(
				functional, minimizer_gsl, tmin);
		FunctionalMinimizer<std::vector<double>, NLopt_vector> fmin_nlopt(
				functional, minimizer_nlopt, tmin);
		switch (MinLib) {
		case gnuscientificlibrary:
			inner_iterations = fmin_gsl(Ndirs, TolFun, tmin);
			break;
		case nonlinearoptimization:
			inner_iterations = fmin_nlopt(Ndirs, TolFun, tmin);
			break;
		default:
			throw MinimizationIllegalValue_exception()
					<< MinimizationIllegalValue_name("MinLib");
			break;
		}
	}
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
		angle_tuple.replace("iteration", (int) (istate.NumberOuterIterations));
		// calculate bregman angles for angles database
		{
			const IterationState::angles_t angles =
					istate.calculateBregmanAngles(newdir);
			for (unsigned int i = 0; (i < MAXANGLES) && (i < angles.size());
					++i) {
				std::stringstream componentname;
				componentname << "bregman_angle" << i + 1;
				angle_tuple.replace(componentname.str(), angles[i]);
			}
		}
		// calculate "scalar product" angles for angles database
		{
			const IterationState::angles_t angles = istate.calculateAngles(
					newdir);
			for (unsigned int i = 0; (i < MAXANGLES) && (i < angles.size());
					++i) {
				std::stringstream componentname;
				componentname << "angle" << i + 1;
				angle_tuple.replace(componentname.str(), angles[i]);
			}
		}
	}
}

void SequentialSubspaceMinimizer::updateSearchspace(
		const SpaceElement_ptr_t& _truesolution,
		const SpaceElement_ptr_t& dual_solution,
		const SpaceElement_ptr_t& newdir,
		const double alpha)
{
	/// update search space with new direction
	if (_truesolution->isZero()) {
		istate.updateSearchSpace(newdir, alpha, dual_solution,
				istate.m_solution);
	} else {
		istate.updateSearchSpace(newdir, alpha, dual_solution, _truesolution);
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
		*dual_x -= tempelement;
	}
	BOOST_LOG_TRIVIAL(trace)<< "x^*_n+1 is " << dual_x;
	*istate.m_solution = refs.J_q(dual_x);
	BOOST_LOG_TRIVIAL(trace)<< "x_n+1 is " << istate.m_solution;
	*_x = istate.m_solution;
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
				residual,
				residuum,
				N);
	}

	/// -# calculate some values prior to loop
	// Jx=DualityMapping(x,NormX,PowerX,TolX);
	SpaceElement_ptr_t dual_solution =
			calculateDualStartingValue(_dualstartvalue);

	// create Bregman distance object
	boost::shared_ptr<BregmanDistance> Delta_p;
	if (!_truesolution->isZero())
		Delta_p.reset(new BregmanDistance (
				refs.NormX,
				dynamic_cast<const PowerTypeDualityMapping &>(refs.J_p),
				refs.J_p.getPower()));

	/// build data tuple for iteration, overall, and angles information
	Table& per_iteration_table = database.addTable("per_iteration");
	Table::Tuple_t per_iteration_tuple = addInfoToPerIterationTable(refs);
	Table& overall_table = database.addTable("overall");
	Table::Tuple_t overall_tuple = addInfoToOverallTable(refs);
	Table& angle_table = database.addTable("angles");
	Table::Tuple_t angle_tuple = addInfoToAnglesTable(refs);

	/// -# check initial stopping criterion
	const double ynorm = refs.NormY(refs.y);
	bool StopCriterion = false;
	const double initial_relative_residuum = fabs(istate.residuum/ynorm);
	StopCriterion = CheckRelativeResiduum(istate.residuum, ynorm);

	while (!StopCriterion) {
		/// Calculation of search direction
		// Jw=DualityMapping(w,NormY,PowerY,TolX);
		searchdir.update(refs, istate.m_residual);
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
		updateSearchspace(_truesolution, dual_solution, searchdir.u, alpha);

		/// get optimal stepwidth
		std::vector<double> tmin(N, 0.);
		double stepwidth_norm = 0.;
		unsigned int inner_iterations = 0;
		try {
			inner_iterations =
					calculateStepWidth(refs, dual_solution, tmin,
							istate.getSearchSpace(), istate.getAlphas());
			stepwidth_norm = std::inner_product(tmin.begin(), tmin.end(), tmin.begin(), stepwidth_norm);
		} catch (MinimizerIllegalNumber_exception &e) {
			BOOST_LOG_TRIVIAL(error)
					<< "Encountered illegal number in line search minimum, not updating.";
			tmin = std::vector<double>(N,0.);
		}

		/// database update prior to iterate update
		per_iteration_tuple.replace( "iteration", (int)istate.NumberOuterIterations);
		per_iteration_tuple.replace( "residual", istate.residuum);
		per_iteration_tuple.replace( "relative_residual", istate.residuum/ynorm);
		per_iteration_tuple.replace( "bregman_distance",
				calculateBregmanDistance(
								Delta_p, istate.m_solution, _truesolution, dual_solution));
		per_iteration_tuple.replace( "error",
				calculateError(istate.m_solution, _truesolution));
		per_iteration_tuple.replace( "updated_index", (int)istate.searchspace->getIndex());
		per_iteration_tuple.replace("inner_iterations",
				(int) (inner_iterations));
		per_iteration_tuple.replace( "stepwidth", sqrt(stepwidth_norm));
		// submit current tuples
		per_iteration_table.addTuple(per_iteration_tuple);
		if (DoCalculateAngles)
			angle_table.addTuple(angle_tuple);

		/// update iterate
		updateIterates(refs, tmin, _problem->x, dual_solution);

		/// update residual
		istate.residuum = calculateResidual(_problem,istate.m_residual);

		/// check iterations count/wall time
		boost::chrono::high_resolution_clock::time_point timing_intermediate =
				boost::chrono::high_resolution_clock::now();
		++istate.NumberOuterIterations;
		StopCriterion =
				CheckIterations(istate.NumberOuterIterations)
				|| CheckRelativeResiduum(istate.residuum, ynorm)
				|| CheckWalltime(boost::chrono::duration<double>(timing_intermediate - timing_start));

		/// check for non-convergence
		const double current_relative_residuum = fabs(istate.residuum/ynorm);
		if (isNonConverging(current_relative_residuum,
				initial_relative_residuum)) {
			fillPerIterationTable(per_iteration_tuple, per_iteration_table);
			fillAngleTable(angle_tuple,angle_table);
		}

		// print intermediate solution
		printIntermediateSolution(
				istate.m_solution, refs.A, istate.NumberOuterIterations);
	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();

	// submit overall_tuple
	overall_tuple.replace( "iterations", istate.NumberOuterIterations );
	overall_tuple.replace( "residual", istate.residuum );
	overall_tuple.replace( "relative_residual", istate.residuum/ynorm );
	overall_tuple.replace( "runtime",
			boost::chrono::duration<double>(timing_end - timing_start).count() );
	finalizeOverallTuple(overall_tuple, refs);
	overall_table.addTuple(overall_tuple);

	// and return solution
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

