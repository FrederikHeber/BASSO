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
#include <boost/bind.hpp>
#include <boost/chrono.hpp>
#include <cassert>
#include <fenv.h>
#include <Minimizations/Mappings/DualityMappingFactory.hpp>
#include <fstream>
#include <sstream>
#include "Log/Logging.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/InverseProblems/QuickAccessReferences.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriteriaFactory.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Norms/Specifics/RegularizedL1Norm.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"
#include "Options/CommandLineOptions.hpp"

//#define BREGMANDISTANCEERRORTHRESHOLD 1

using namespace boost::assign;

// static entities
GeneralMinimizer::MinLib_names_t GeneralMinimizer::MinLib_names;

GeneralMinimizer::GeneralMinimizer(
		const CommandLineOptions &_opts,
		const InverseProblem_ptr_t &_inverseproblem,
		Database &_database
		) :
	dbcontainer(_database),
	Delta(_opts.delta),
	MaxOuterIterations(_opts.maxiter),
	MaxInnerIterations(_opts.maxinneriter),
	TolX(_opts.tolerance_spacex),
	TolY(Delta),
	TolFun(_opts.tolerance_linesearch),
	outputsteps(_opts.outputsteps),
	MinLib(gnuscientificlibrary),
	OldBregmanDistance(0.),
	l2norm(NormFactory::getInstance().create(
			"lp",
			_inverseproblem->x->getSpace(),
			NormFactory::args_t(1, boost::any(2.))))
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

	// create stopping criterion
	StoppingCriteriaFactory stop_factory;
	const_cast<StoppingCriterion::ptr_t &>(stopping_criteria) =
			stop_factory.create(_opts.stopping_criteria, _opts.stopping_args);

	// set NoOp additional parameters
	// don't do this in initializer list as class not fully constructed
	dbcontainer.setAddParamsCallback(
			boost::bind(&GeneralMinimizer::addAdditionalParametersToTuple,
					boost::cref(*this), _1, _2));
}

GeneralMinimizer::~GeneralMinimizer()
{
	// add view if not present to database if not empty
	if (dbcontainer.size() != 0)
		if (!dbcontainer.createViews())
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
		u = _refs.A_t(Jw);
		if (u->getSpace()->getDimension() > 10)
			BOOST_LOG_TRIVIAL(trace)
					<< "newdir is " << u;
		else
			BOOST_LOG_TRIVIAL(debug)
					<< "newdir is " << u;
}

bool GeneralMinimizer::CheckStoppingCondition(
		const boost::chrono::duration<double> &_time,
		const int _current_outeriterations,
		const double _residuum,
		const double _ynorm) const
{
	const bool result =
			(*stopping_criteria)(
					_time, _current_outeriterations, _residuum, _ynorm);
	return result;
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
	<< "dual_x_n is " << m_dual_solution;
	BOOST_LOG_TRIVIAL(trace)
	<< "R_n is " << m_residual;
}

double GeneralMinimizer::calculateResidual(
		const InverseProblem_ptr_t &_problem,
		SpaceElement_ptr_t &_residual
		) const
{
	const Mapping &A = static_cast<const Mapping &>(*_problem->A);
	*_residual = A(_problem->x);
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
		const Mapping &_A,
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
					ost << _A(_solution);
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Could not write all data of projected intermediate solution to stream.\n";
				}
			else {
				std::cerr << "Failed to open " << solution_file.str() << std::endl;
			}
		}

	}
}

void GeneralMinimizer::resetState()
{
	// reset Bregman distance check variable
	resetBregmanDistance();
	// reset state of derived objects
	resetState_interal();
}
