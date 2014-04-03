/*
 * LandweberMinimizer.cpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#include "LandweberMinimizer.hpp"

#include <algorithm>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <levmar.h>
#include <limits>

#include "DualityMapping.hpp"
#include "LpNorm.hpp"
#include "MinimizationExceptions.hpp"
#include "SmoothnessModulus.hpp"

LandweberMinimizer::LandweberMinimizer() :
	MaxOuterIterations(50),
	TolX(1e-6),
	TolFun(1e-12),
	C(0.5)
{}

struct LandweberParameters
{
	SmoothnessModulus *modul;
	double lambda;
};

/** Static function to wrap call to BregmanFunctional::operator()().
 *
 */
static void
func(double *p, double *hx, int m, int n, void *adata)
{
	LandweberParameters *params = static_cast<LandweberParameters *>(adata);
	const double result = (*params->modul)(p[0]);
	hx[0] = result/p[0] - params->lambda;
}

double LandweberMinimizer::calculateMatchingTau(
		SmoothnessModulus &_modul,
		const double _lambda
		) const
{
	LandweberParameters params;
	params.modul = &_modul;
	params.lambda = _lambda;
	double tau = .5;
	double info[LM_INFO_SZ];
	double x[1];
	const int m = 1;
	const int n = 1;
	double lb[] = { std::numeric_limits<double>::epsilon() };
	double ub[] = { 1. };
	double *work = (double *)malloc( ( LM_DIF_WORKSZ(m, n) + m*m) * sizeof(double));
	double *covar = work+LM_DIF_WORKSZ(m, n);
	int ret = dlevmar_bc_dif(
			*func,
			&tau,
			x,
			1,
			1,
			lb,
			ub,
			NULL,
			1000,
			NULL, 	/* opts[4] */
			info,
			work,
			covar,
			&params
			);
	if (ret == -1)
		throw MinimizationFunctionError_exception()
			<< MinimizationFunctionError_name(tau);
	// free everything
	free(work);
	BOOST_LOG_TRIVIAL(trace)
		<< "Matching tau from modulus of smoothness is " << tau;
	BOOST_LOG_TRIVIAL(trace)
		<< "Counter-check: rho(tau)/tau = "
		<< _modul(tau)/tau << ", compare with lambda = " << _lambda;

	return tau;
}

GeneralMinimizer::ReturnValues
LandweberMinimizer::operator()(
		const Eigen::VectorXd &_x0,
		const unsigned int _val_NormX,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y,
		const unsigned int _val_NormY,
		const double _PowerY,
		const double _Delta
		) const
{
	double PowerX;
	double DualPowerX;
	double val_DualNormX = _val_NormX/(_val_NormX - 1.);
	double G;
	if (_val_NormX > 2) {
		PowerX = _val_NormX;
		DualPowerX = PowerX/(PowerX - 1.);
		G = ::pow(2., 2. - val_DualNormX);
	} else {
		PowerX = 2.;
		DualPowerX = 2.;
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

	BOOST_LOG_TRIVIAL(debug) << "Calculating "
			<< _A << "*" << _x0 << "-" << _y;
	Eigen::VectorXd w = _A * _x0 - _y;
	double wNorm = NormY(w);
	double TolY = _Delta;
	returnvalues.residuum = wNorm;

	// check stopping criterion
	StopCriterion = (fabs(returnvalues.residuum) <= TolY);

	// start building up search space 'U' with the search vectors 'u'
	DualityMapping J_p(_val_NormX);
	DualityMapping j_r(_val_NormY);
	J_p.setTolerance(TolX);
	j_r.setTolerance(TolY);
	DualityMapping J_q(val_DualNormX);
	J_q.setTolerance(TolX);
	Eigen::VectorXd dual_solution =
			J_p(returnvalues.solution, PowerX);
	while (!StopCriterion) {
		BOOST_LOG_TRIVIAL(debug)
				<< "Current iteration is " << returnvalues.NumberOuterIterations
				<< " at position " << returnvalues.solution.transpose();
				DualityMapping J_y(_val_NormY);
		J_y.setTolerance(TolX);
		Eigen::VectorXd u = _A.transpose()*J_y(w, _PowerY);
		double alpha = 0.;
		const double _ANorm = _A.norm();
		if ((returnvalues.NumberOuterIterations == 0)
			|| (dual_solution.isZero())) {
			alpha =
					C * ::pow(_val_NormY, _val_NormX - 1.) / ::pow(_ANorm,_val_NormX)
					* ::pow(returnvalues.residuum, _val_NormX - _val_NormY);
		} else {
			SmoothnessModulus modul(_val_NormX);
			const double modulus = modul(1);
			alpha = C/(::pow(2., val_DualNormX) * G * _ANorm)
					* returnvalues.residuum / dual_solution.norm();
			const double lambda = std::min(modulus, alpha);
			// find intermediate value in smoothness modulus to match lambda
			const double tau = calculateMatchingTau(modul, lambda);
			// calculate step width
			alpha = (tau/_ANorm)
					* ::pow( dual_solution.norm(), _val_NormX-1.)
					/ ::pow(returnvalues.residuum, _val_NormY - 1.);
		}

		// iterate
		dual_solution -= alpha * _A.transpose()
				* j_r( _A * returnvalues.solution - _y, _PowerY);

		// finally map back from X^{\conj} to X
		returnvalues.solution =
				J_q(dual_solution, PowerX);

		// check iterations count
		++returnvalues.NumberOuterIterations;
		StopCriterion =
				(returnvalues.NumberOuterIterations >= MaxOuterIterations);
	}

	return returnvalues;
}
