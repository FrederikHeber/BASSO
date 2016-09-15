/*
 * SequentialSubspaceMinimizerNoise.cpp
 *
 *  Created on: Jun 03, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "SequentialSubspaceMinimizerNoise.hpp"

#include <algorithm>
#include <boost/bind.hpp>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <iterator>
#include <numeric>
#include <sstream>

#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "Log/Logging.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/InverseProblems/QuickAccessReferences.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/BregmanProjectionFunctional.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer_deprecated.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Functions/Minimizers/Minimizer.hpp"
#include "Minimizations/Functions/Minimizers/MinimizerExceptions.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Norms/Norm.hpp"

SequentialSubspaceMinimizerNoise::SequentialSubspaceMinimizerNoise(
		const CommandLineOptions &_opts,
		const InverseProblem_ptr_t &_inverseproblem,
		Database &_database
		) :
	SequentialSubspaceMinimizer(
		_opts,
		_inverseproblem,
		_database
		)
{}

SequentialSubspaceMinimizerNoise::ReturnValues
SequentialSubspaceMinimizerNoise::operator()(
		const InverseProblem_ptr_t &_problem,
		const SpaceElement_ptr_t &_startvalue,
		const SpaceElement_ptr_t &_dualstartvalue,
		const SpaceElement_ptr_t &_truesolution
		)
{
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	// make sure we have two search directions
	// TODO: Remove this when scheme allows for any number of search directions
	assert ( N <= 2);

	// gather some refs for convenient access
	QuickAccessReferences refs(_problem);

	// G constant
	const double G =
			refs.NormX.getPvalue() > 2. ?
			::pow(2., 2. - refs.DualNormX.getPvalue()) :
			 refs.DualNormX.getPvalue() - 1.;
	LOG(debug, "G is " << G);

	/// initialize return structure
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
	dbcontainer.setParameterKey(
			refs.NormX.getPvalue(),
			refs.NormY.getPvalue(),
			N,
			refs.SpaceX.getDimension(),
			MaxOuterIterations);
	Table::Tuple_t& per_iteration_tuple = dbcontainer.preparePerIterationTuple();
	per_iteration_tuple.insert( std::make_pair("inner_iterations", (int)0), Table::Data);
	Table::Tuple_t& overall_tuple = dbcontainer.prepareOverallTuple();

	/// -# check stopping criterion
	const double ynorm = refs.NormY(refs.y);
	const double initial_residuum = istate.residuum;
	bool StopCriterion = CheckStoppingCondition(
			boost::chrono::duration<double>(0.),
			istate.NumberOuterIterations,
			istate.residuum,
			ynorm);

	std::vector<double> d(N,0.);
	istate.status = ReturnValues::started;
	unsigned int tuple_counter = 1;
	double stepwidth_norm = 0.;
	unsigned int inner_iterations = 0;
	while (!StopCriterion) {
		/// Calculation of search direction
		searchdir.update(refs, istate.m_residual);

		/// output prior to iterate update
		istate.output(ynorm);

		// uNorm=norm(u,DualNormX);
		const double uNorm = refs.DualNormX(searchdir.u);
		LOG(debug, "uNorm is " << uNorm);
		const double ScaleFactor = 1./uNorm;
		// alpha=u'*x-Residual^PowerY; (29) or (24) in original form
		const double alpha =
				searchdir.u * istate.m_solution - ::pow(istate.residuum, refs.j_r.getPower());
		LOG(trace, "alpha is " << alpha);
		updateSearchspace(_truesolution, searchdir.u, alpha);
		const unsigned int index = istate.searchspace->getIndex();

		// d=Delta*Residual^(PowerY-1); (24) or (32)
		d[index] = Delta * ::pow(istate.residuum, refs.j_r.getPower()-1.);
		// beta=Residual^(PowerY-1)*(Residual-Delta)/uNorm^DualPowerX;
		const double beta =
				::pow(istate.residuum, refs.j_r.getPower()-1.)
				 *(istate.residuum-Delta)/::pow(uNorm, refs.J_q.getPower());

		/// database update prior to iterate update
		updatePerIterationTuple(per_iteration_tuple,
				ynorm, stepwidth_norm, inner_iterations, Delta_p, _truesolution);

		/// 1st update onto new hyperplane
		// index contains newest direction on whose associated hyperplane
		// x in general does not lie.
		try {
			std::vector<double> tmin(1, 0.);
			if (istate.m_dual_solution->isApproxToConstant(0, TolX)) {
				// tmin=beta^(PowerX-1);
				tmin[0] = ::pow(beta, refs.J_p.getPower() - 1.);
				LOG(trace, "tmin is " << tmin[0]);
			} else {
				// t0=(beta/G)^(PowerX-1);
				// wrong, rather we have
				// t0=(Residual*(Residual-Delta)/(G*uNorm^DualPowerX))^(PowerX-1);
				// see Remark b) on page 17
				std::vector<double> t0(1, 0.);
				tmin[0] = uNorm * ::pow(beta/G,refs.J_p.getPower() - 1.);
				LOG(trace, "Initial tmin[0] is " << t0[0]);

				std::vector<SpaceElement_ptr_t> ConstrainedSearchSpace(1);
				ConstrainedSearchSpace[0] =
						ScaleFactor * istate.getSearchSpace()[index];
				std::vector<double> steps(1);
				steps[0] = ScaleFactor*(istate.getAlphas()[index]+d[index]);

				inner_iterations =
						calculateStepWidth(refs, istate.m_dual_solution, tmin,
								ConstrainedSearchSpace, steps);

				tmin[0] *= ScaleFactor;
			}
			stepwidth_norm = std::inner_product(tmin.begin(), tmin.end(), tmin.begin(), stepwidth_norm);
			// x=DualityMapping(Jx-tmin*u,DualNormX,DualPowerX,TolX);
			/// update iterate
			tmin.resize(2, 0.);
			if (index != 0)
				std::swap(tmin[index], tmin[0]);
			updateIterates(refs, tmin, _problem->x, istate.m_dual_solution);
		} catch (MinimizerIllegalNumber_exception &e) {
			LOG(error, "Encountered illegal number in line search minimum, not updating.");
		}

		if ((istate.NumberOuterIterations > 0) && (N == 2)) {
			/// 2nd update onto intersection of hyperspaces
			// find "old" index (one before current index)
			const std::vector<unsigned int> &lastIndices = istate.searchspace->getLastIndices();
			const std::vector<unsigned int>::const_iterator iter =
					std::find( lastIndices.begin(), lastIndices.end(), (unsigned int)1);
			// with N == 2 we must always
			assert ( iter != lastIndices.end() );
			const unsigned int lastIndex =
					std::distance( lastIndices.begin(), iter );
			const double Resold = istate.getSearchSpace()[lastIndex] * istate.m_solution;

			// project onto correct second hyperplane intersection
			std::vector<double> steps(N);
			steps[index] = istate.getAlphas()[index]+d[index];
			steps[lastIndex] = istate.getAlphas()[lastIndex];
			LOG(debug, "On intersection? "
				<< Resold << " in ["
				<< steps[lastIndex]-d[lastIndex]
				<< "," << steps[lastIndex]+d[lastIndex]<< "]?");
			// numerically stabler: first subtract hyperplane offset
			const double planeOffset = Resold - steps[lastIndex];
			if (fabs(planeOffset) > d[lastIndex]) {
				if (planeOffset > d[lastIndex]) {
					steps[lastIndex] = steps[lastIndex]+d[lastIndex];
				} else {
					steps[lastIndex] = steps[lastIndex]-d[lastIndex];
				}

				// rescale noisy offset
				std::transform(
						steps.begin(), steps.end(),
						steps.begin(),
						std::bind2nd(std::multiplies<double>(), ScaleFactor));


				// scale search space such that current dir has norm 1
				std::vector<SpaceElement_ptr_t> ScaledSearchSpace(N, SpaceElement_ptr_t());
				for (size_t i=0;i<N;++i)
					ScaledSearchSpace[i] = ScaleFactor * istate.getSearchSpace()[i];

				std::vector<double> tmin(N, 0.);
//				const unsigned int inner_iterations =
						calculateStepWidth(refs, istate.m_dual_solution, tmin,
								ScaledSearchSpace, steps);

				// scale tmin back
				std::transform(
						tmin.begin(), tmin.end(),
						tmin.begin(),
						std::bind2nd(std::multiplies<double>(), ScaleFactor));

				stepwidth_norm = std::inner_product(tmin.begin(), tmin.end(), tmin.begin(), stepwidth_norm);

				/// update iterate
				updateIterates(refs, tmin, _problem->x, istate.m_dual_solution);
			}
		}
		per_iteration_tuple.replace( "stepwidth", sqrt(stepwidth_norm));

		/// submit current tuples
		if (everynthtuple != 0) {
			if (tuple_counter >= everynthtuple) {
				dbcontainer.data_per_iteration_table.addTuple(per_iteration_tuple);
				tuple_counter = 1;
			} else
				++tuple_counter;
		}

		/// update residual
		istate.residuum = calculateResidual(_problem,istate.m_residual);

		/// check iterations count/wall time
		boost::chrono::high_resolution_clock::time_point timing_intermediate =
				boost::chrono::high_resolution_clock::now();
		++istate.NumberOuterIterations;
		const double current_residuum = istate.residuum;
		StopCriterion = CheckStoppingCondition(
			timing_intermediate - timing_start,
			istate.NumberOuterIterations,
			istate.residuum,
			ynorm);

		/// check for non-convergence
		if (isNonConverging(current_residuum,
				initial_residuum))
			fillPerIterationTable(per_iteration_tuple, tuple_counter);

		/// print intermediate solution
		printIntermediateSolution(
				istate.m_solution, refs.A, istate.NumberOuterIterations);
	}
	/// last output of istate
	istate.output(ynorm);

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	LOG(debug, "Iteration was stopped by "
			<< stopping_criteria->whoIsTrue(
					timing_end - timing_start,
					istate.NumberOuterIterations,
					istate.residuum,
					ynorm));

	// submit last per_iteration tuple
	if (everynthtuple != 0) {
		updatePerIterationTuple(per_iteration_tuple,
				ynorm, stepwidth_norm, inner_iterations, Delta_p, _truesolution);
		dbcontainer.data_per_iteration_table.addTuple(per_iteration_tuple);
	}

	// submit overall_tuple
	overall_tuple.replace( "iterations", istate.NumberOuterIterations );
	overall_tuple.replace( "residual", istate.residuum );
	if (fabs(ynorm) > BASSOTOLERANCE)
		overall_tuple.replace( "relative_residual", istate.residuum/ynorm );
	else
		overall_tuple.replace( "relative_residual", 0. );
	overall_tuple.replace( "runtime",
			boost::chrono::duration<double>(timing_end - timing_start).count() );
	dbcontainer.finalizeOverallTuple(overall_tuple, refs);
	dbcontainer.data_overall_table.addTuple(overall_tuple);

	// and return solution
	istate.status = ReturnValues::finished;
	return istate;
}
