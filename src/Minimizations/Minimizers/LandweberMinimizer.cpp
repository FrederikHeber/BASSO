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
 * LandweberMinimizer.cpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LandweberMinimizer.hpp"

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/chrono.hpp>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>

#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/InverseProblems/QuickAccessReferences.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/ResidualFunctional.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/RelativeShrinkageMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/StepWidths/DetermineStepWidth.hpp"
#include "Minimizations/Minimizers/StepWidths/DetermineStepWidthFactory.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/Specifics/RegularizedL1Norm.hpp"
#include "Options/CommandLineOptions.hpp"

LandweberMinimizer::LandweberMinimizer(
		const CommandLineOptions &_opts,
		const InverseProblem_ptr_t &_inverseproblem,
		Database &_database
		) :
	GeneralMinimizer(
			_opts,
			_inverseproblem,
			_database
			),
	C(0.9),
	stepwidth_type((const enum DetermineStepWidthFactory::stepwidth_enumeration)_opts.stepwidth_type)
{}

void LandweberMinimizer::setC(const double _C)
{
	if ((_C <= 0.)) // || ( _C > 1.))
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("C");

	const_cast<double&>(C) = _C;
}

void LandweberMinimizer::updatePerIterationTuple(
		Table::Tuple_t& per_iteration_tuple,
		ReturnValues &returnvalues,
		const double &ynorm,
		const double &alpha,
		boost::shared_ptr<BregmanDistance> &Delta_p,
		const SpaceElement_ptr_t &_truesolution) const
{
	per_iteration_tuple.replace( "iteration", (int)returnvalues.NumberOuterIterations);
	per_iteration_tuple.replace( "residual", returnvalues.residuum);
	if (fabs(ynorm) > BASSOTOLERANCE)
		per_iteration_tuple.replace( "relative_residual", returnvalues.residuum/ynorm);
	else
		per_iteration_tuple.replace( "relative_residual", 0.);
	per_iteration_tuple.replace( "bregman_distance",
			calculateBregmanDistance(
					Delta_p, returnvalues.m_solution, _truesolution, returnvalues.m_dual_solution));
	per_iteration_tuple.replace( "error",
			calculateError(returnvalues.m_solution, _truesolution));
	per_iteration_tuple.replace( "stepwidth", alpha);
}

GeneralMinimizer::ReturnValues
LandweberMinimizer::operator()(
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
	ReturnValues returnvalues;
	returnvalues.status = ReturnValues::starting;
	returnvalues.NumberOuterIterations = 0;
	// set iterate 'x' as start vector 'x0'
	returnvalues.m_solution = refs.SpaceX.createElement();
	*returnvalues.m_solution = _startvalue;
	returnvalues.m_dual_solution = refs.DualSpaceX.createElement();
	*returnvalues.m_dual_solution = _dualstartvalue;
	returnvalues.m_residual = refs.SpaceY.createElement();
	// calculate starting residual and norm
	returnvalues.residuum = calculateResidual(
			_problem,
			returnvalues.m_residual);
	const double ynorm = refs.NormY(refs.y);

	/// set initial lambda
	double mutual_coherence = 0.;
	bool setLambdaAdaptively = false;
	const RegularizedL1Norm * const regl1_norm = dynamic_cast<const RegularizedL1Norm *>(&refs.NormX);
	if ((stepwidth_type >=
			DetermineStepWidthFactory::ConstantRegularizedL1Norm)
			&& (regl1_norm != NULL)
			&& (regl1_norm->getLambda() == 0.)) {
		mutual_coherence = dynamic_cast<const LinearMapping &>(refs.A).MutualCoherence();
		setRegularizationParameter(mutual_coherence, _truesolution);
		setLambdaAdaptively = true;
	}

	/// build data tuple for iteration and overall information
	dbcontainer.setParameterKey(
			refs.NormX.getPvalue(),
			refs.NormY.getPvalue(),
			1,
			refs.SpaceX.getDimension(),
			MaxOuterIterations);
	Table::Tuple_t& per_iteration_tuple = dbcontainer.preparePerIterationTuple();
	Table::Tuple_t& overall_tuple = dbcontainer.prepareOverallTuple();
//	overall_tuple.insert( std::make_pair("runtime_matrix_vector_products", (int)0), Table::Data );
//	overall_tuple.insert( std::make_pair("runtime_vector_vector_products", (int)0), Table::Data );
	// due to Eigen's lazy evaluation runtime is not measured accurately

	/// -# check stopping criterion
	bool StopCriterion = false;
	StopCriterion = CheckStoppingCondition(
			boost::chrono::duration<double>(0.),
			returnvalues.NumberOuterIterations,
			returnvalues.residuum,
			ynorm);

	// create Bregman distance object
	boost::shared_ptr<BregmanDistance> Delta_p;
	if (!_truesolution->isZero())
		Delta_p.reset(new BregmanDistance (
				refs.NormX, refs.J_p, refs.J_p.getPower()));

	// step width algorithm
	ResidualFunctional::calculateResidual_t residualizer =
			boost::bind(&LandweberMinimizer::calculateResidual,
					boost::cref(*this), _1, _2);
	DetermineStepWidth_ptr_t stepwidth =
			DetermineStepWidthFactory::createInstance(
					_problem,
					stepwidth_type,
					C,
					residualizer,
					_problem->y->getSpace()->getDualityMapping());

	/// -# loop over stopping criterion
	returnvalues.status = ReturnValues::started;
	unsigned int tuple_counter = 1;
	double alpha = 0.;
	while (!StopCriterion) {
		/// Update adjoint (if non-linear
		refs.A.updateAdjoint(returnvalues.m_solution);

		/// Calculation of search direction
		searchdir.update(refs, returnvalues.m_residual);

		/// output prior to iterate update
		returnvalues.output(ynorm);

		/// find step width
		// (F. Schöpfer, 11.4.2014) too conservative! Line search instead
		alpha = (*stepwidth)(
				returnvalues.m_dual_solution,
				searchdir.u,
				returnvalues.m_solution,
				returnvalues.m_residual,
				returnvalues.residuum,
				TolX,
				0.
				);
		LOG(trace, "Step width is " << alpha);

		/// database update prior to iterate update
		updatePerIterationTuple(per_iteration_tuple,
				returnvalues, ynorm, alpha, Delta_p, _truesolution);

		/// update iterate
		*returnvalues.m_dual_solution -= alpha * searchdir.u;
		LOG(trace, "x^*_n+1 is " << returnvalues.m_dual_solution);
		// finally map back from X^{\conj} to X: x_{n+1}
		*returnvalues.m_solution = refs.J_q(returnvalues.m_dual_solution);
		LOG(trace, "x_n+1 is " << returnvalues.m_solution);
		*_problem->x = returnvalues.m_solution;
		if (setLambdaAdaptively) {
			setRegularizationParameter(
					mutual_coherence,
					_truesolution);
		}

		/// update residual
		returnvalues.residuum = calculateResidual(
				_problem,
				returnvalues.m_residual);

		/// check iterations count/wall time
		boost::chrono::high_resolution_clock::time_point timing_intermediate =
				boost::chrono::high_resolution_clock::now();
		++returnvalues.NumberOuterIterations;
		StopCriterion = CheckStoppingCondition(
			timing_intermediate - timing_start,
			returnvalues.NumberOuterIterations,
			returnvalues.residuum,
			ynorm);

		/// print intermediate solution
		printIntermediateSolution(
				returnvalues.m_solution, refs.A, returnvalues.NumberOuterIterations);

		/// submit current tuple to database
		if (everynthtuple != 0) {
			if (tuple_counter >= everynthtuple) {
				dbcontainer.data_per_iteration_table.addTuple(per_iteration_tuple);
				tuple_counter = 1;
			} else
				++tuple_counter;
		}
	}
	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	LOG(debug, "Iteration was stopped by "
			<< stopping_criteria->whoIsTrue(
					timing_end - timing_start,
					returnvalues.NumberOuterIterations,
					returnvalues.residuum,
					ynorm));

	// submit last per_iteration tuple
	if (everynthtuple != 0) {
		updatePerIterationTuple(per_iteration_tuple,
				returnvalues, ynorm, alpha, Delta_p, _truesolution);
		dbcontainer.data_per_iteration_table.addTuple(per_iteration_tuple);
	}

	// submit overall_tuple
	overall_tuple.replace( "iterations", returnvalues.NumberOuterIterations );
	overall_tuple.replace( "residual", returnvalues.residuum );
	if (fabs(ynorm) > BASSOTOLERANCE)
		overall_tuple.replace( "relative_residual", returnvalues.residuum/ynorm );
	else
		overall_tuple.replace( "relative_residual", 0. );
	overall_tuple.replace( "runtime",
			boost::chrono::duration<double>(timing_end - timing_start).count() );
	dbcontainer.finalizeOverallTuple(overall_tuple, refs);
	dbcontainer.data_overall_table.addTuple(overall_tuple);

	returnvalues.status = ReturnValues::finished;

	return returnvalues;
}

void LandweberMinimizer::setRegularizationParameter(
		const double mutual_coherence,
		const SpaceElement_ptr_t &_solution) const
{
	// count the number of zeros in _solution
	unsigned int nonzerocomponents = 0;
	for (unsigned int i=0;i<_solution->getSpace()->getDimension();++i)
		if (fabs((*_solution)[i]) > std::numeric_limits<double>::epsilon())
			++nonzerocomponents;
	// strict positivity criterion
	if ( nonzerocomponents < .5*(1. + 1./mutual_coherence) ) {
		const double lambda =
				1./(1. + 1./mutual_coherence - 2. * (double)nonzerocomponents);
		const RegularizedL1Norm &regularizednorm =
				dynamic_cast<const RegularizedL1Norm &>(
						*_solution->getSpace()->getNorm()
						);
		const_cast<RegularizedL1Norm &>(regularizednorm).setLambda(lambda);
		LOG(trace, "Lambda of NormX is now " << regularizednorm.getLambda());
		const RelativeShrinkageMapping &mapping =
				dynamic_cast<const RelativeShrinkageMapping &>(
						*_solution->getSpace()->getDualSpace()->getDualityMapping()
						);
		const_cast<RelativeShrinkageMapping &>(mapping).setLambda(lambda);
		LOG(trace, "Lambda of RelativeShrinkageMapping in X^* is now " << mapping.getLambda());
	} else {
		LOG(error, "Cannot set lambda adaptively as criterion is not fulfilled.");
	}
}
