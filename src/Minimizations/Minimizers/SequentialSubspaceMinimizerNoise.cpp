/*
 * SequentialSubspaceMinimizerNoise.cpp
 *
 *  Created on: Jun 03, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "SequentialSubspaceMinimizerNoise.hpp"

#include <boost/log/trivial.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <iterator>
#include <numeric>
#include <sstream>

#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/InverseProblems/QuickAccessReferences.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/BregmanProjectionFunctional.hpp"
#include "Minimizations/Functions/HyperplaneProjection.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Functions/Minimizers/Minimizer.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Norms/Norm.hpp"

SequentialSubspaceMinimizerNoise::SequentialSubspaceMinimizerNoise(
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		const unsigned int _maxinneriter,
		Database &_database,
		const unsigned int _outputsteps
		) :
	SequentialSubspaceMinimizer(
		_inverseproblem,
		_Delta,
		_maxiter,
		_maxinneriter,
		_database,
		_outputsteps
		),
	tau(1.1)
{}

void SequentialSubspaceMinimizerNoise::setTau(
		const double _tau
		)
{
	// check that regularization parameter is greater than 1
	if (_tau <= 1.)
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("tau");
	const_cast<double&>(tau) = _tau;
	// change y tolerance according to regularization parameter
	const_cast<double&>(TolY) = tau * Delta;
}

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

	// gather some refs for convenient access
	QuickAccessReferences refs(_problem);

	// G constant used in theoretical step width
	const double G =
			refs.NormX.getPvalue() > 2. ?
			::pow(2., 2. - refs.DualNormX.getPvalue()) :
			 refs.DualNormX.getPvalue() - 1.;
	BOOST_LOG_TRIVIAL(debug)
		<< "G is " << G;

	// reset inner state if problem has changed
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

	// calculate some values prior to loop
	// Jx=DualityMapping(x,NormX,PowerX,TolX);
	SpaceElement_ptr_t dual_solution = refs.DualSpaceX.createElement();
	*dual_solution = _dualstartvalue;
	BOOST_LOG_TRIVIAL(trace)
		<< "Jx_0 is " << dual_solution;

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

	/// -# check stopping criterion
	const double ynorm = refs.NormY(refs.y);
	bool StopCriterion = false;
	const double initial_residuum = istate.residuum;
	StopCriterion = CheckResiduum(initial_residuum);

	while (!StopCriterion) {
		per_iteration_tuple.replace( "iteration", (int)istate.NumberOuterIterations);
		if ((istate.NumberOuterIterations == 0)
			|| (istate.residuum > TolY)) {
			/// Calculation of search direction
			searchdir.update(refs, istate.m_residual);

			/// output prior to iterate update
			istate.output(ynorm);

			/// update search space with new direction
			// uNorm=norm(u,DualNormX);
			const double uNorm = refs.DualNormX(searchdir.u);
			BOOST_LOG_TRIVIAL(debug)
				<< "uNorm is " << uNorm;
			// alpha=u'*x-Residual^PowerY;
			const double alpha =
					searchdir.u * istate.m_solution - ::pow(istate.residuum,refs.j_r.getPower());
			BOOST_LOG_TRIVIAL(trace)
				<< "alpha is " << alpha;
			updateSearchspace(_truesolution, dual_solution, searchdir.u, alpha);

			// d=Delta*Residual^(PowerY-1);
			const double d =
					Delta * ::pow(istate.residuum,(double)refs.j_r.getPower()-1.);
			// beta=Residual^(PowerY-1)*(Residual-Delta)/uNorm^DualPowerX;
			const double beta =
					::pow(istate.residuum,(double)refs.j_r.getPower()-1.)
					* (istate.residuum-Delta)/::pow(uNorm, refs.J_q.getPower());

			/// database update prior to iterate update
			per_iteration_tuple.replace( "iteration", (int)istate.NumberOuterIterations);
			per_iteration_tuple.replace( "relative_residual", istate.residuum/ynorm);
			per_iteration_tuple.replace( "bregman_distance",
					calculateBregmanDistance(
							Delta_p, istate.m_solution, _truesolution, dual_solution));
			per_iteration_tuple.replace( "error",
					calculateError(istate.m_solution, _truesolution));
			per_iteration_tuple.replace( "updated_index", (int)istate.searchspace->getIndex());
//			per_iteration_tuple.replace("inner_iterations",
//					(int) (inner_iterations));
//			per_iteration_tuple.replace( "stepwidth", sqrt(stepwidth_norm));

			std::vector<double> tmin(1, 0.);
			if (dual_solution->isApproxToConstant(0, TolX)) {
				// tmin=beta^(PowerX-1);
				tmin[0] = ::pow(beta, refs.J_p.getPower() - 1.);
				BOOST_LOG_TRIVIAL(debug)
					<< "tmin is " << tmin[0];
			} else {
				// t0=(beta/G)^(PowerX-1);
				std::vector<double> t0(1, 0.);
				t0[0] = ::pow(beta/G, refs.J_p.getPower() - 1.);
				BOOST_LOG_TRIVIAL(trace)
					<< "t0[0] is " << t0[0];
				{
					// tmin=fminunc(@(t) BregmanProjectionFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
					BregmanProjectionFunctional bregman(
							refs.DualNormX,
							dynamic_cast<const PowerTypeDualityMapping &>(refs.J_q),
							refs.J_q.getPower());

					const HyperplaneProjection functional(
							bregman,
							dual_solution,
							istate.getSearchSpace(),
							istate.getAlphas());
					Minimizer<gsl_vector> minimizer(1);

					// TODO: current alpha needs to be modified for minimization!
					std::vector<double> steps(1);
					steps[0] = alpha+d;
					FunctionalMinimizer<std::vector<double>, gsl_vector> fmin(
						functional, minimizer, t0);

					tmin = t0;
					const unsigned int iterations =
							fmin(1, TolY, tmin);

					std::stringstream output_stepwidth;
					std::copy(tmin.begin(), tmin.end(), std::ostream_iterator<double>(output_stepwidth, " "));
					BOOST_LOG_TRIVIAL(debug)
						<< "tmin " << output_stepwidth.str()
						<< " found in " << iterations << " iterations.";
				}
			}
			double stepwidth_norm = 0.;
			stepwidth_norm = std::inner_product(tmin.begin(), tmin.end(), tmin.begin(), stepwidth_norm);
			per_iteration_tuple.replace( "stepwidth", sqrt(stepwidth_norm));
			// x=DualityMapping(Jx-tmin*u,DualNormX,DualPowerX,TolX);
			/// update iterate
			updateIterates(refs, tmin, _problem->x, dual_solution);

			/// update residual
			istate.residuum = calculateResidual(
					_problem,
					istate.m_residual);

			/// check iterations count/wall time
			boost::chrono::high_resolution_clock::time_point timing_intermediate =
					boost::chrono::high_resolution_clock::now();
			++istate.NumberOuterIterations;
			const double current_residuum = istate.residuum;
			StopCriterion =
					CheckIterations(istate.NumberOuterIterations)
					|| CheckResiduum(current_residuum)
					|| CheckWalltime(boost::chrono::duration<double>(timing_intermediate - timing_start));

			/// check for non-convergence
			if (isNonConverging(current_residuum,
					initial_residuum))
				fillPerIterationTable(per_iteration_tuple, per_iteration_table);

			/// print intermediate solution
			printIntermediateSolution(
					istate.m_solution, refs.A, istate.NumberOuterIterations);
		}


		/// submit current tuple
		per_iteration_table.addTuple(per_iteration_tuple);
	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	std::cout << "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start)
			<< "." << std::endl;

	// submit overall_tuple
	overall_tuple.replace( "iterations", istate.NumberOuterIterations );
	overall_tuple.replace( "relative_residual", istate.residuum );
	overall_tuple.replace( "runtime",
			boost::chrono::duration<double>(timing_end - timing_start).count() );
	finalizeOverallTuple(overall_tuple, refs);
	overall_table.addTuple(overall_tuple);

	// and return solution
	return istate;
}
