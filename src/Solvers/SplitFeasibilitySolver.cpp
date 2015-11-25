/*
 * SplitFeasibilitySolver.cpp
 *
 *  Created on: Nov 19, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SplitFeasibilitySolver.hpp"

#include <cassert>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "Log/Logging.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Options/CommandLineOptions.hpp"
#include "Solvers/AuxiliaryConstraintsProblem.hpp"

SplitFeasibilitySolver::SplitFeasibilitySolver(
		const CommandLineOptions &_opts
		) :
		opts(_opts),
		name("SplitFeasibilityProblem")
{}

void SplitFeasibilitySolver::registerFeasibilityProblem(
		FeasibilityProblem::ptr_t &_fp)
{
	problems.push_back(_fp);
}

void SplitFeasibilitySolver::registerAuxiliaryConstraints(
		const AuxiliaryConstraints::ptr_t &_auxiliary_constraints
		)
{
	FeasibilityProblem::ptr_t fp(
			new AuxiliaryConstraintsProblem(_auxiliary_constraints));
	registerFeasibilityProblem(fp);
}

static std::string createBar(const unsigned int _length)
{
	std::vector<std::string> barvector(_length, "=");
	std::stringstream output;
	std::copy(barvector.begin(), barvector.end(),
			std::ostream_iterator<std::string>(output, ""));
	return output.str();
}

GeneralMinimizer::ReturnValues SplitFeasibilitySolver::operator()(
		const SpaceElement_ptr_t &_startingvalue)
{
	GeneralMinimizer::ReturnValues result;
	result.m_solution = _startingvalue->getSpace()->createElement();
	*result.m_solution = _startingvalue;
	result.m_dual_solution =
			_startingvalue->getSpace()->getDualSpace()->createElement();
	result.m_dual_solution->setZero();
	result.residuum = std::numeric_limits<double>::max();
	result.status = GeneralMinimizer::ReturnValues::notbegun;

	BOOST_LOG_TRIVIAL(info)
			<< "Starting SplitFeasibilityProblem ...";

	const std::string bar = createBar(30);
	for (unsigned int SplitFeasibilityProblem_loops = 0;
			(result.residuum > opts.delta)
			&& (SplitFeasibilityProblem_loops < opts.max_sfp_loops);
			++SplitFeasibilityProblem_loops) {
		BOOST_LOG_TRIVIAL(info)
				<< bar << " SFP #" << SplitFeasibilityProblem_loops
				<< " " << bar;

		// go through each problem
		for (problems_t::iterator iter = problems.begin();
				iter != problems.end();
				++iter) {
			(*iter)->clear();
			result = (**iter)(result.m_solution);
			(*iter)->finish();
			BOOST_LOG_TRIVIAL(info)
					<< "Residual after problem " << (*iter)->getName()
					<< " is " << result.residuum;
			if (result.status != GeneralMinimizer::ReturnValues::finished) {
				BOOST_LOG_TRIVIAL(error)
						<< "The last SFP part did not finish properly.";
				break;
			}
		}
		BOOST_LOG_TRIVIAL(info)
				<< bar
				<< createBar(6+ceil(SplitFeasibilityProblem_loops/10))
				<< bar;
	}

	BOOST_LOG_TRIVIAL(info)
			<< "Finishing SplitFeasibilityProblem ...";

	return result;
}
