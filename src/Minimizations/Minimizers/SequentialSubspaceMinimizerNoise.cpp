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

#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Functions/BregmanProjectionFunctional.hpp"
#include "Minimizations/Functions/FunctionMinimizer.hpp"
#include "Minimizations/Functions/HyperplaneProjection.hpp"
#include "Minimizations/Functions/MinimizationFunctional.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Norms/Norm.hpp"

// instantiate required template functions
CONSTRUCT_FUNCTIONMINIMIZER(Eigen::VectorXd)

SequentialSubspaceMinimizerNoise::SequentialSubspaceMinimizerNoise(
		const InverseProblem_ptr_t &_inverseproblem,
		const double _Delta,
		const unsigned int _maxiter,
		Database &_database,
		const unsigned int _outputsteps
		) :
	SequentialSubspaceMinimizer(
		_inverseproblem,
		_Delta,
		_maxiter,
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
	const NormedSpace & SpaceX = *_problem->A->getSourceSpace();
	const NormedSpace & DualSpaceX = *SpaceX.getDualSpace();
	const NormedSpace & SpaceY = *_problem->A->getTargetSpace();
	const Norm & NormX = *SpaceX.getNorm();
	const Norm & DualNormX = *DualSpaceX.getNorm();
	const Norm & NormY = *SpaceY.getNorm();
	const Mapping & J_p = *SpaceX.getDualityMapping();
	const Mapping & J_q = *DualSpaceX.getDualityMapping();
	const Mapping & j_r = *SpaceY.getDualityMapping();
	const SpaceElement_ptr_t &y = _problem->y;
	const LinearMapping &A =
			dynamic_cast<const LinearMapping &>(*_problem->A);
	const Mapping_ptr_t A_adjoint = A.getAdjointMapping();
	const LinearMapping &A_t =
			dynamic_cast<const LinearMapping &>(*A_adjoint);

	// G constant used in theoretical step width
	const double G =
			NormX.getPvalue() < 2. ?
			::pow(2., 2. - DualNormX.getPvalue()) :
			 DualNormX.getPvalue() - 1.;
	BOOST_LOG_TRIVIAL(trace)
		<< "G is " << G;

	/// -# initialize return structure
	ReturnValues returnvalues;
	returnvalues.NumberOuterIterations = 0;
	// set iterate 'x' as start vector 'x0'
	returnvalues.m_solution = _problem->x->getSpace()->createElement();
	*returnvalues.m_solution = _startvalue;
	returnvalues.m_residual = _problem->y->getSpace()->createElement();
// calculate starting residual and norm
	returnvalues.residuum = calculateResidual(
			_problem,
			returnvalues.m_residual);
	const double ynorm = NormY(y);

	/// -# check stopping criterion
	bool StopCriterion = false;
	StopCriterion = (fabs(returnvalues.residuum/ynorm) <= TolY);

	// calculate some values prior to loop
	// Jx=DualityMapping(x,NormX,PowerX,TolX);
	SpaceElement_ptr_t dual_solution = DualSpaceX.createElement();
	*dual_solution = _dualstartvalue;
	BOOST_LOG_TRIVIAL(trace)
		<< "Jx_0 is " << dual_solution;

	// build data tuple for iteration information
	Table& per_iteration_table = database.addTable("per_iteration");
	Table::Tuple_t per_iteration_tuple;
	per_iteration_tuple.insert( std::make_pair("p", NormX.getPvalue()));
	per_iteration_tuple.insert( std::make_pair("r", NormY.getPvalue()));
	per_iteration_tuple.insert( std::make_pair("dim", (int)SpaceX.getDimension()));
	per_iteration_tuple.insert( std::make_pair("iteration", (int)0));
	per_iteration_tuple.insert( std::make_pair("stepwidth", (int)0));
	per_iteration_tuple.insert( std::make_pair("relative_residual", 0.));
	per_iteration_tuple.insert( std::make_pair("error", 0.));
	per_iteration_tuple.insert( std::make_pair("bregman_distance", 0.));

	// build data tuple for overall information
	Table& overall_table = database.addTable("overall");
	Table::Tuple_t overall_tuple;
	overall_tuple.insert( std::make_pair("p", NormX.getPvalue()));
	overall_tuple.insert( std::make_pair("r", NormY.getPvalue()));
	overall_tuple.insert( std::make_pair("dim", (int)SpaceX.getDimension()));
	overall_tuple.insert( std::make_pair("iterations", (int)0));
	overall_tuple.insert( std::make_pair("relative_residual", 0.));
	overall_tuple.insert( std::make_pair("runtime", 0.));
	overall_tuple.insert( std::make_pair("matrix_vector_products", (int)0));
	overall_tuple.insert( std::make_pair("runtime_matrix_vector_products", 0.));
	overall_tuple.insert( std::make_pair("vector_vector_products", (int)0));
	overall_tuple.insert( std::make_pair("runtime_vector_vector_products", 0.));

	// reset inner state of problem has changed
	if (state.getDimension() != A.getSourceSpace()->getDimension())
		state.set(A.getSourceSpace()->getDimension(), N);
	while (!StopCriterion) {
		per_iteration_tuple.replace( "iteration", (int)returnvalues.NumberOuterIterations);
		if ((returnvalues.NumberOuterIterations == 0)
			|| (returnvalues.residuum > TolY)) {
			BOOST_LOG_TRIVIAL(debug)
					<< "#" << returnvalues.NumberOuterIterations
					<< " with residual of " << returnvalues.residuum;
			BOOST_LOG_TRIVIAL(trace)
					<< "x_n is " << returnvalues.m_solution;
			BOOST_LOG_TRIVIAL(trace)
					<< "R_n is " << returnvalues.m_residual;
			per_iteration_tuple.replace( "relative_residual", returnvalues.residuum/ynorm);

			const SpaceElement_ptr_t Jw = j_r( returnvalues.m_residual );
			BOOST_LOG_TRIVIAL(trace)
					<< "j_r (residual) is " << Jw;

			// u=A'*DualityMapping(w,NormY,PowerY,TolX);
			const SpaceElement_ptr_t u = A_t * Jw;
			BOOST_LOG_TRIVIAL(trace)
				<< "u is " << u;

			// uNorm=norm(u,DualNormX);
			const double uNorm = DualNormX(u);
			BOOST_LOG_TRIVIAL(trace)
				<< "uNorm is " << uNorm;
			// alpha=u'*x-Residual^PowerY;
			const double alpha =
					u * returnvalues.m_solution - ::pow(returnvalues.residuum,j_r.getPower());
			BOOST_LOG_TRIVIAL(trace)
				<< "alpha is " << alpha;
			// d=Delta*Residual^(PowerY-1);
			const double d =
					Delta * ::pow(returnvalues.residuum,(double)j_r.getPower()-1.);
			// beta=Residual^(PowerY-1)*(Residual-Delta)/uNorm^DualPowerX;
			const double beta =
					::pow(returnvalues.residuum,(double)j_r.getPower()-1.)
					* (returnvalues.residuum-Delta)/::pow(uNorm, J_q.getPower());

			Eigen::VectorXd tmin(1);
			tmin.setZero();
			if (dual_solution->isApproxToConstant(0, TolX)) {
				// tmin=beta^(PowerX-1);
				tmin[0] = ::pow(beta, J_p.getPower() - 1.);
				BOOST_LOG_TRIVIAL(trace)
					<< "tmin is " << tmin;
			} else {
				// t0=(beta/G)^(PowerX-1);
				Eigen::VectorXd t0(1);
				t0[0] = ::pow(beta/G, J_p.getPower() - 1.);
				BOOST_LOG_TRIVIAL(trace)
					<< "t0[0] is " << t0[0];
				{
					// tmin=fminunc(@(t) BregmanProjectionFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
					BregmanProjectionFunctional bregman(
							DualNormX,
							dynamic_cast<const PowerTypeDualityMapping &>(J_q),
							J_q.getPower(),
							MatrixVectorProduct_subspace,
							ScalarVectorProduct_subspace);

					HyperplaneProjection functional(
							bregman,
							dual_solution->getVectorRepresentation(),
							state.U,
							state.alphas);

					// TODO: current alpha needs to be modified for minimization!
					Eigen::VectorXd steps(1);
					steps[0] = alpha+d;
					FunctionMinimizer<Eigen::VectorXd> minimizer(
						functional, t0);

					tmin = minimizer(1, TolY, tmin);

					BOOST_LOG_TRIVIAL(trace)
						<< "tmin is " << tmin.transpose();
				}
			}
			per_iteration_tuple.replace( "stepwidth", tmin.norm());
//			const Eigen::VectorXd uold = u;
//			const double alphao0ld = alpha;
//			const double dold = d;
			// x=DualityMapping(Jx-tmin*u,DualNormX,DualPowerX,TolX);
			{
				const SpaceElement_ptr_t tempelement = DualSpaceX.createElement();
				*tempelement = tmin*u->getVectorRepresentation();
				*dual_solution -= tempelement;
			}
			BOOST_LOG_TRIVIAL(trace)
					<< "x^*_n+1 is " << dual_solution;
			returnvalues.m_solution = J_q(dual_solution);
			BOOST_LOG_TRIVIAL(trace)
					<< "x_n+1 is " << returnvalues.m_solution;
			*_problem->x = returnvalues.m_solution;

			// update residual
			returnvalues.residuum = calculateResidual(
					_problem,
					returnvalues.m_residual);

			// check iterations count
			++returnvalues.NumberOuterIterations;
			StopCriterion =
					(returnvalues.NumberOuterIterations >= MaxOuterIterations)
					|| (fabs(returnvalues.residuum) <= TolY);

			// print intermediat solution
			printIntermediateSolution(
					 returnvalues.m_solution->getVectorRepresentation(),
					 A.getMatrixRepresentation(),
					returnvalues.NumberOuterIterations);
		}

		// submit current tuple
		per_iteration_table.addTuple(per_iteration_tuple);
	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	std::cout << "The operation took " << boost::chrono::duration<double>(timing_end - timing_start) << "." << std::endl;

	// submit overall_tuple
	overall_tuple.replace( "iterations", returnvalues.NumberOuterIterations );
	overall_tuple.replace( "relative_residual", returnvalues.residuum );
	overall_tuple.replace( "runtime",
			boost::chrono::duration_cast<boost::chrono::duration<double> >(timing_end - timing_start).count() );
	overall_tuple.replace( "matrix_vector_products", (int)MatrixVectorProduct.getCount() );
	overall_tuple.replace( "vector_vector_products", (int)ScalarVectorProduct.getCount() );
	overall_tuple.replace( "matrix_vector_products_subspace", (int)MatrixVectorProduct_subspace.getCount() );
	overall_tuple.replace( "vector_vector_products_subspace", (int)ScalarVectorProduct_subspace.getCount() );
	overall_table.addTuple(overall_tuple);

	// and return solution
	return returnvalues;
}
