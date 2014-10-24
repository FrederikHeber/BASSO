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

#include "Minimizations/Mappings/LpDualityMapping.hpp"
#include "Norms/LpNorm.hpp"
#include "MinimizationExceptions.hpp"
#include "Minimizations/Functions/BregmanProjectionFunctional.hpp"
#include "Minimizations/Functions/FunctionMinimizer.hpp"
#include "Minimizations/Functions/HyperplaneProjection.hpp"
#include "Minimizations/Functions/MinimizationFunctional.hpp"

// instantiate required template functions
CONSTRUCT_FUNCTIONMINIMIZER(Eigen::VectorXd)

SequentialSubspaceMinimizerNoise::SequentialSubspaceMinimizerNoise(
		const DualityMappingsContainer &_container,
		const double _NormY,
		const double _PowerY,
		const double _Delta,
		const unsigned int _maxiter,
		Database &_database,
		const unsigned int _outputsteps
		) :
	SequentialSubspaceMinimizer(
		_container,
		_NormY,
		_PowerY,
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
		const Eigen::VectorXd &_x0,
		const Eigen::VectorXd &_dualx0,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y,
		const Eigen::VectorXd &_solution
		)
{
//	NoCols = _A.innerSize();
//	NoRows = _A.outerSize();

	double PowerX;
	double DualPowerX;
	double G;
	if (val_NormX > 2) {
		PowerX = val_NormX;
		DualPowerX = PowerX/(PowerX - 1.);
		G = ::pow(2., 2. - val_DualNormX);
	} else {
		PowerX = 2.;
		DualPowerX = 2.;
		G = val_DualNormX - 1.;
	}
	BOOST_LOG_TRIVIAL(trace)
		<< "New PowerX is " << PowerX;
	BOOST_LOG_TRIVIAL(trace)
		<< "New DualPowerX is " << DualPowerX;
	BOOST_LOG_TRIVIAL(trace)
		<< "G is " << G;

	/// -# initialize return structure
	ReturnValues returnvalues;
	returnvalues.NumberOuterIterations = 0;
	// set iterate 'x' as start vector 'x0'
	returnvalues.solution = _x0;
	// calculate starting residual and norm
	returnvalues.residuum = calculateResidual(
			_x0, _A, _y,
			returnvalues.residual);

	/// -# check stopping criterion
	bool StopCriterion = false;
	StopCriterion = (fabs(returnvalues.residuum) <= TolY);

	// calculate some values prior to loop
	// Jx=DualityMapping(x,NormX,PowerX,TolX);
	Eigen::VectorXd dual_solution = _dualx0;
	BOOST_LOG_TRIVIAL(trace)
		<< "Jx_0 is " << dual_solution.transpose();
//	const double modulus_at_one = modul(1);
//	const double _ANorm = ::pow(2, 1.+ 1./val_NormY); //_A.norm();
//	BOOST_LOG_TRIVIAL(trace)
//		<< "_ANorm " << _ANorm;

	// reset inner state of problem has changed
	if (state.getDimension() != _A.outerSize())
		state.set(_A.outerSize(), N);
	while (!StopCriterion) {
		if ((returnvalues.NumberOuterIterations == 0)
			|| (returnvalues.residuum > TolY)) {
			BOOST_LOG_TRIVIAL(debug)
					<< "#" << returnvalues.NumberOuterIterations
					<< " with residual of " << returnvalues.residuum;
			BOOST_LOG_TRIVIAL(trace)
					<< "x_n is " << returnvalues.solution.transpose();
			BOOST_LOG_TRIVIAL(trace)
					<< "R_n is " << returnvalues.residual.transpose();
			BOOST_LOG_TRIVIAL(trace)
				<< "j_r (residual) is "
				<< j_r( returnvalues.residual, PowerY).transpose();

			// u=A'*DualityMapping(w,NormY,PowerY,TolX);
			Eigen::VectorXd u =
					_A.transpose()*j_r(returnvalues.residual, PowerY);
			BOOST_LOG_TRIVIAL(trace)
				<< "u is " << u.transpose();

			// uNorm=norm(u,DualNormX);
			const double uNorm = DualNormX(u);
			BOOST_LOG_TRIVIAL(trace)
				<< "uNorm is " << uNorm;
			// alpha=u'*x-Residual^PowerY;
			const double alpha =
					u.transpose() * _x0 - ::pow(returnvalues.residuum,PowerY);
			BOOST_LOG_TRIVIAL(trace)
				<< "alpha is " << alpha;
			// d=Delta*Residual^(PowerY-1);
			const double d =
					Delta * ::pow(returnvalues.residuum,(double)PowerY-1.);
			// beta=Residual^(PowerY-1)*(Residual-Delta)/uNorm^DualPowerX;
			const double beta =
					::pow(returnvalues.residuum,(double)PowerY-1.)
					* (returnvalues.residuum-Delta)/::pow(uNorm, DualPowerX);

			Eigen::VectorXd tmin(1);
			tmin.setZero();
			if (dual_solution.isApproxToConstant(0, TolX)) {
				// tmin=beta^(PowerX-1);
				tmin[0] = ::pow(beta, PowerX - 1.);
				BOOST_LOG_TRIVIAL(trace)
					<< "tmin is " << tmin;
			} else {
				// t0=(beta/G)^(PowerX-1);
				Eigen::VectorXd t0(1);
				t0[0] = ::pow(beta/G, PowerX - 1.);
				BOOST_LOG_TRIVIAL(trace)
					<< "t0[0] is " << t0[0];
				{
					// tmin=fminunc(@(t) BregmanProjectionFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
					BregmanProjectionFunctional bregman(
							DualNormX,
							J_q,
							MatrixVectorProduct_subspace,
							ScalarVectorProduct_subspace);

					HyperplaneProjection functional(
							bregman,
							dual_solution,
							state.U,
							state.alphas,
							DualPowerX);

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
//			const Eigen::VectorXd uold = u;
//			const double alphao0ld = alpha;
//			const double dold = d;
			// x=DualityMapping(Jx-tmin*u,DualNormX,DualPowerX,TolX);
			dual_solution -= tmin*u;
			BOOST_LOG_TRIVIAL(trace)
					<< "x^*_n+1 is " << dual_solution.transpose();
			returnvalues.solution =
					J_q(dual_solution , DualPowerX);
			BOOST_LOG_TRIVIAL(trace)
					<< "x_n+1 is " << returnvalues.solution.transpose();

			// update residual
			returnvalues.residuum = calculateResidual(
					returnvalues.solution,
					_A,
					_y,
					returnvalues.residual);

			// check iterations count
			++returnvalues.NumberOuterIterations;
			StopCriterion =
					(returnvalues.NumberOuterIterations >= MaxOuterIterations)
					|| (fabs(returnvalues.residuum) <= TolY);

			// print intermediat solution
			printIntermediateSolution(
					returnvalues.solution,
					_A,
					returnvalues.NumberOuterIterations);
		}
	}

	// and return solution
	return returnvalues;
}
