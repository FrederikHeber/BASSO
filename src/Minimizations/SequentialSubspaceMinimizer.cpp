/*
 * SequentialSubspaceMinimizer.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "SequentialSubspaceMinimizer.hpp"

#include <boost/log/trivial.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <levmar.h>

#include "BregmanFunctional.hpp"
#include "DualityMapping.hpp"
#include "LpNorm.hpp"
#include "MinimizationExceptions.hpp"

SequentialSubspaceMinimizer::SequentialSubspaceMinimizer(
		const unsigned int _NormX,
		const unsigned int _NormY,
		const double _PowerX,
		const double _PowerY,
		const double _Delta,
		const unsigned int _maxiter,
		const unsigned int _outputsteps
		) :
	GeneralMinimizer(
			_NormX,
			_NormY,
			_PowerX,
			_PowerY,
			_Delta,
			_maxiter,
			_outputsteps
			),
	tau(1.1)
{}

void SequentialSubspaceMinimizer::setTau(
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

/** Structure containing all parameters to call BregmanFunctional functions.
 *
 * This is required to use levmar function minimization that only allows
 * to pass a void* pointer to pass on information to the function to be
 * minimized.
 *
 * \sa BregmanFunctional
 *
 */
struct BregmanParameters
{
	BregmanFunctional &bregman;
	Eigen::VectorXd t;
	const Eigen::VectorXd &x;
	const Eigen::MatrixXd &U;
	const Eigen::VectorXd &alpha;
	const double q;

	/** Constructor to initialize refs.
	 *
	 */
	BregmanParameters(
		BregmanFunctional &_bregman,
		const Eigen::VectorXd &_t,
		const Eigen::VectorXd &_x,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha,
		const double _q
		) :
			bregman(_bregman),
			t(_t),
			x(_x),
			U(_U),
			alpha(_alpha),
			q(_q)
	{}
};

/** Static function to wrap call to BregmanFunctional::operator()().
 *
 */
static void
func(double *p, double *hx, int m, int n, void *adata)
{
	struct BregmanParameters *params =
			static_cast<BregmanParameters *>(adata);
	params->t[0] = p[0];
	hx[0] =
			(params->bregman)(
					params->t,
					params->x,
					params->U,
					params->alpha,
					params->q);
}

/** Static function to wrap call to BregmanFunctional::operator()().
 *
 */
static void
jacf(double *p, double *j, int m, int n, void *adata)
{
	struct BregmanParameters *params =
			static_cast<BregmanParameters *>(adata);
	// update where to evaluate
	params->t[0] = p[0];
	Eigen::VectorXd grad =
			(params->bregman).gradient(
					params->t,
					params->x,
					params->U,
					params->alpha,
					params->q);
	for (int i=0; i<m;++i)
		j[i] = grad[i];
}

SequentialSubspaceMinimizer::ReturnValues
SequentialSubspaceMinimizer::operator()(
		const Eigen::VectorXd &_x0,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y
		) const
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
	Eigen::VectorXd dual_solution =
			J_p(returnvalues.solution, PowerX);
	BOOST_LOG_TRIVIAL(trace)
		<< "Jx_0 is " << dual_solution.transpose();
//	const double modulus_at_one = modul(1);
//	const double _ANorm = ::pow(2, 1.+ 1./val_NormY); //_A.norm();
//	BOOST_LOG_TRIVIAL(trace)
//		<< "_ANorm " << _ANorm;

	// start building up search space 'U' with the search vectors 'u'
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

			double tmin = 0.;
			if (dual_solution.isApproxToConstant(0, TolX)) {
				// tmin=beta^(PowerX-1);
				tmin = ::pow(beta, PowerX - 1.);
				BOOST_LOG_TRIVIAL(trace)
					<< "tmin is " << tmin;
			} else {
				// t0=(beta/G)^(PowerX-1);
				Eigen::VectorXd t0(1);
				t0[0] = ::pow(beta/G, PowerX - 1.);
				BOOST_LOG_TRIVIAL(trace)
					<< "t0[0] is " << t0[0];
				// tmin=fminunc(@(t) BregmanFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
				BregmanFunctional bregman(val_DualNormX, TolX);
				// TODO_ we actually have to perform a function minimization
				// with respect to t starting at t0
				Eigen::VectorXd steps(1);
				steps[0] = alpha+d;
				BregmanParameters params(
						bregman,
						t0,
						dual_solution,
						u,
						steps,
						DualPowerX);
				double info[LM_INFO_SZ];
				double x[1];
				const int m = 1;
				const int n = 1;
				double *work = (double *)malloc( ( LM_DIF_WORKSZ(m, n) + m*m) * sizeof(double));
				double *covar = work+LM_DIF_WORKSZ(m, n);
				int ret = dlevmar_der(
						*func,
						*jacf,
						&tmin,
						x,
						m,
						n,
						1000, 	/* itmax */
						NULL, 	/* opts[4] */
						info,
						work,
						covar,
						static_cast<void *>(&params)
						);
				if (ret == -1)
					throw MinimizationFunctionError_exception()
						<< MinimizationFunctionError_name(tau);
				// free everything
				free(work);
				// solution is in tmin already
				BOOST_LOG_TRIVIAL(trace)
					<< "tmin is " << tmin;
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
