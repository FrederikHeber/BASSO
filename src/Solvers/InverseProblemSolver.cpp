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
 * InverseProblemSolver.cpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InverseProblemSolver.hpp"

#include <cassert>

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixFactorizer/Helpers/detail.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/SingularValueDecomposition.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriteriaFactory.hpp"
#include "Options/CommandLineOptions.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"

InverseProblemSolver::InverseProblemSolver(
		InverseProblem_ptr_t &_inverseproblem,
		Database_ptr_t &_database,
		const CommandLineOptions &_opts,
		const bool _checkTrueSolution
		) :
		inverseproblem(_inverseproblem),
		database(_database),
		opts(_opts),
		checkTrueSolution(_checkTrueSolution),
		name("InverseProblem")
{
	minimizer = SolverFactory::createMinimizer(
					opts, _inverseproblem, database);
}

GeneralMinimizer::ReturnValues InverseProblemSolver::operator()(
		const SpaceElement_ptr_t &_startingvalue,
		const SpaceElement_ptr_t _truesolution)
{
	GeneralMinimizer::ReturnValues result;
	result.residuum = 0.;

	if (minimizer == NULL) {
		LOG(error, "Minimizer could not be constructed, exiting.");
		result.status = GeneralMinimizer::ReturnValues::error;
		return result;
	}

	SpaceElement_ptr_t truesolution;
	if (checkTrueSolution) {
		if (_truesolution) {
			truesolution = _truesolution;
		} else {
			const LinearMapping &A = static_cast<LinearMapping&>(*inverseproblem->A);
			const SpaceElement_ptr_t &rhs = inverseproblem->y;
			SingularValueDecomposition svd = A.getSVD();
			const SpaceElement_ptr_t truesolution = svd.solve(rhs);

			// empty or true solution from diagonalization
			LOG(trace, "True solution is " << *truesolution
					<< " with norm " << (A(truesolution) - rhs)->Norm()/rhs->Norm());
		}
	} else {
		truesolution =
				inverseproblem->x->getSpace()->createElement();
	}

	// prepare start value and dual solution
	result.m_solution = _startingvalue->getSpace()->createElement();
	*result.m_solution = _startingvalue;
	*inverseproblem->x = result.m_solution;
	if (result.m_solution->getSpace()->getDimension() < 10) {
		LOG(debug, "Starting at x0 = " << result.m_solution);
	}

	// only for smooth spaces we may use the duality mapping
	if (inverseproblem->x->getSpace()->getNorm()->isSmooth()) {
		result.m_dual_solution =
				(*inverseproblem->x->getSpace()->getDualityMapping())(result.m_solution);
	} else {
		result.m_dual_solution =
				inverseproblem->x->getSpace()->getDualSpace()->createElement();
	}

	/// solver inverse problem
	try{
		result = (*minimizer)(
						inverseproblem,
						result.m_solution,
						result.m_dual_solution,
						truesolution);
		minimizer->resetState();
	} catch (MinimizationIllegalValue_exception &e) {
		std::cerr << "Illegal value for "
				<< *boost::get_error_info<MinimizationIllegalValue_name>(e)
				<< std::endl;
		result.status = GeneralMinimizer::ReturnValues::error;
		return result;
	}
	assert( *result.m_solution == *inverseproblem->x );
	assert( result.m_solution->getSpace().get() == inverseproblem->x->getSpace().get() );

	return result;
}

void InverseProblemSolver::clear()
{
	database->clear();
}

void InverseProblemSolver::finish()
{
	database->finish();
}
