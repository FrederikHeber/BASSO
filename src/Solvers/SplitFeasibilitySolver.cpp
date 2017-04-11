/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
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
		const SpaceElement_ptr_t &_startingvalue,
		const SpaceElement_ptr_t)
{
	GeneralMinimizer::ReturnValues result;
	result.m_solution = _startingvalue->getSpace()->createElement();
	*result.m_solution = _startingvalue;
	result.m_dual_solution =
			_startingvalue->getSpace()->getDualSpace()->createElement();
	result.residuum = std::numeric_limits<double>::max();
	result.status = GeneralMinimizer::ReturnValues::notbegun;

	LOG(debug, "Starting SplitFeasibilityProblem ...");

	const std::string bar = createBar(30);
	for (unsigned int SplitFeasibilityProblem_loops = 0;
			(result.residuum > opts.delta)
			&& (SplitFeasibilityProblem_loops < opts.max_sfp_loops);
			++SplitFeasibilityProblem_loops) {
		LOG(debug, bar << " SFP #" << SplitFeasibilityProblem_loops << " " << bar);

		// go through each problem
		for (problems_t::iterator iter = problems.begin();
				iter != problems.end();
				++iter) {
			(*iter)->clear();
			result = (**iter)(result.m_solution);
			(*iter)->finish();
			LOG(debug, "Residual after problem " << (*iter)->getName() << " is " << result.residuum);
			if (result.status != GeneralMinimizer::ReturnValues::finished) {
				LOG(error, "The last SFP part did not finish properly.");
				break;
			}
		}
		LOG(debug, bar
				<< createBar(6+ceil(SplitFeasibilityProblem_loops/10)) << bar);
	}

	LOG(debug, "Finishing SplitFeasibilityProblem ...");

	return result;
}
