/*
 * SequentialSubspaceMinimizer.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "SequentialSubspaceMinimizer.hpp"

#include <cassert>
#include <cmath>
#include <Eigen/Dense>
#include <levmar.h>

#include "BregmanFunctional.hpp"
#include "DualityMapping.hpp"
#include "LpNorm.hpp"

SequentialSubspaceMinimizer::SequentialSubspaceMinimizer() :
	MaxOuterIterations(1000),
	TolX(1e-6),
	TolFun(1e-12),
	tau(1.1)
{
	// check that regularization parameter is greater than 1
	assert( tau > 1 );
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
		const unsigned int _val_NormX,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y,
		const unsigned int _val_NormY,
		const double _PowerY,
		const double _Delta
		)
{
//	NoCols = _A.innerSize();
//	NoRows = _A.outerSize();

	double PowerX;
	double DualPowerX;
	double val_DualNormX;
	double G;
	if (_val_NormX > 2) {
		PowerX = _val_NormX;
		DualPowerX = PowerX/(PowerX - 1.);
		G = ::pow(2., 2. - val_DualNormX);
	} else {
		PowerX = 2.;
		DualPowerX = 2.;
		val_DualNormX = _val_NormX/(_val_NormX - 1.);
		G = val_DualNormX - 1.;
	}

	// prepare norm functors
	LpNorm NormX(_val_NormX);
	LpNorm NormY(_val_NormY);
	LpNorm DualNormX(val_DualNormX);

	// initialize return structure
	ReturnValues returnvalues;
	returnvalues.NumberOuterIterations = 0;
	// set iterate 'x' as start vector 'x0'
	returnvalues.solution = _x0;

	bool StopCriterion = false;

	Eigen::VectorXd w = _A * _x0 - _y;
	double wNorm = NormY(w);
	double TolY = tau * _Delta;
	returnvalues.residuum = wNorm;

	// check stopping criterion
	StopCriterion = (fabs(returnvalues.residuum) <= TolY);

	// start building up search space 'U' with the search vectors 'u'
	while (!StopCriterion) {
		if ((returnvalues.NumberOuterIterations == 0)
			|| (returnvalues.residuum > TolY)) {
			// u=A'*DualityMapping(w,NormY,PowerY,TolX);
			DualityMapping J_y(_val_NormY);
			J_y.setTolerance(TolX);
			Eigen::VectorXd u = _A.transpose()*J_y(w, _PowerY);
			// uNorm=norm(u,DualNormX);
			const double uNorm = DualNormX(u);
			// alpha=u'*x-Residual^PowerY;
			const double alpha =
					u.transpose() * _x0 - ::pow(returnvalues.residuum,_PowerY);
			// d=Delta*Residual^(PowerY-1);
			const double d = _Delta * ::pow(returnvalues.residuum,(double)_PowerY-1.);
			// beta=Residual^(PowerY-1)*(Residual-Delta)/uNorm^DualPowerX;
			const double beta = ::pow(returnvalues.residuum,(double)_PowerY-1.)
					* (returnvalues.residuum-_Delta)/::pow(uNorm, DualPowerX);
			// Jx=DualityMapping(x,NormX,PowerX,TolX);
			DualityMapping J_x(_val_NormX);
			J_x.setTolerance(TolX);
			Eigen::VectorXd Jx = J_x(_x0, PowerX);

			double tmin;
			if (Jx.isZero(BASSOTOLERANCE)) {
				// tmin=beta^(PowerX-1);
				tmin = ::pow(beta, PowerX - 1.);
			} else {
				// t0=(beta/G)^(PowerX-1);
				Eigen::VectorXd t0(1);
				t0[0] = ::pow(beta/G, PowerX - 1.);
				// tmin=fminunc(@(t) BregmanFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
				BregmanFunctional bregman(val_DualNormX, TolX);
				// TODO_ we actually have to perform a function minimization
				// with respect to t starting at t0
				Eigen::VectorXd steps(1);
				steps[0] = alpha+d;
				BregmanParameters params(
						bregman,
						t0,
						Jx,
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
				// free everything
				free(work);
				free(covar);
				// solution is in tmin already
			}
			const Eigen::VectorXd uold = u;
			const double alphao0ld = alpha;
			const double dold = d;
			// x=DualityMapping(Jx-tmin*u,DualNormX,DualPowerX,TolX);
			DualityMapping Jdual_x(val_DualNormX);
			Jdual_x.setTolerance(TolX);
			returnvalues.solution = Jdual_x(Jx - tmin*u, DualPowerX);

			// check iterations count
			++returnvalues.NumberOuterIterations;
			StopCriterion =
					(returnvalues.NumberOuterIterations >= MaxOuterIterations);
		}
	}

	// and return solution
	return returnvalues;
}