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

LandweberMinimizer::LandweberMinimizer(
		const double _NormX,
		const double _NormY,
		const double _PowerY,
		const double _Delta,
		const unsigned int _maxiter
		) :
	val_NormX(_NormX),
	val_NormY(_NormY),
	PowerY(_PowerY),
	Delta(_Delta),
	MaxOuterIterations(_maxiter),
	TolX(1e-6),
	TolFun(1e-12),
	C(0.9)
{
	if ((C <= 0.)) // || ( C > 1.))
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("C");
}

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
	const double norm = result/p[0] - params->lambda;
	hx[0] = norm*norm;
}

double LandweberMinimizer::calculateMatchingTau(
		SmoothnessModulus &_modul,
		const double _lambda
		) const
{
	LandweberParameters params;
	params.modul = &_modul;
	params.lambda = _lambda;
	double tau[] = { .5 };
	double info[LM_INFO_SZ];
	double x[] = { 0. };
	const int m = 1;
	const int n = 1;
	double lb[] = { std::numeric_limits<double>::epsilon() };
	double ub[] = { 1. };
	double *work = (double *)malloc( ( LM_DIF_WORKSZ(m, n) + m*m) * sizeof(double));
	double *covar = work+LM_DIF_WORKSZ(m, n);
	int ret = dlevmar_bc_dif(
			(*func),
			tau,
			x,
			m,
			n,
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
//	  std::cout << "Minimization info: ";
//	  for(int i=0; i<LM_INFO_SZ; ++i)
//	    std::cout << info[i] << ", ";
//	  std::cout << std::endl;
	  if (ret == -1)
		throw MinimizationFunctionError_exception()
			<< MinimizationFunctionError_name(tau[0]);
	// free everything
	free(work);
	BOOST_LOG_TRIVIAL(trace)
		<< "Matching tau from modulus of smoothness is " << tau[0];
	BOOST_LOG_TRIVIAL(trace)
		<< "Counter-check: rho(tau)/tau = "
		<< _modul(tau[0])/tau[0] << ", compare with lambda = " << _lambda;

	return tau[0];
}

double LandweberMinimizer::calculateResidual(
		const Eigen::VectorXd &_x0,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y,
		const LpNorm &_NormY,
		Eigen::VectorXd &_residual
		) const
{
	_residual = _A * _x0 - _y;
	return _NormY(_residual);
}

GeneralMinimizer::ReturnValues
LandweberMinimizer::operator()(
		const Eigen::VectorXd &_x0,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y
		) const
{
	double PowerX;
	double DualPowerX;
	double val_DualNormX = val_NormX/(val_NormX - 1.);
	double G;
	const double TolY = Delta;
	PowerX = val_NormX;
	DualPowerX = val_DualNormX; // PowerX/(PowerX - 1.);
	if (val_DualNormX < 2.) {
		G = ::pow(2., 2. - val_DualNormX);
	} else {
//		PowerX = 2.;
//		DualPowerX = 2.;
		G = val_DualNormX - 1.;
	}
	BOOST_LOG_TRIVIAL(debug)
		<< "p is " << val_NormX
		<< ", q is " << val_DualNormX
		<< ", r is " << val_NormY
		<< ", power of J_p is " <<  PowerX
		<< ", power of J_q is " <<  DualPowerX
		<< ", power of J_r is " <<  PowerY;

	// prepare norm and duality mapping functors
	LpNorm NormX(val_NormX);
	LpNorm NormY(val_NormY);
	LpNorm DualNormX(val_DualNormX);
	DualityMapping J_p(val_NormX);
	DualityMapping J_q(val_DualNormX);
	DualityMapping j_r(val_NormY);
	J_p.setTolerance(TolX);
	J_q.setTolerance(TolX);
	j_r.setTolerance(TolY);
	SmoothnessModulus modul(val_DualNormX);

	// initialize return structure
	ReturnValues returnvalues;
	returnvalues.NumberOuterIterations = 0;
	// set iterate 'x' as start vector 'x0'
	returnvalues.solution = _x0;

	bool StopCriterion = false;

	BOOST_LOG_TRIVIAL(debug) << "Calculating "
			<< _A << "*" << _x0 << "-" << _y;

	// calculate starting residual and norm
	returnvalues.residuum = calculateResidual(
			_x0, _A, _y, NormY,
			returnvalues.residual);

	// check stopping criterion
	StopCriterion = (fabs(returnvalues.residuum) <= TolY);

	Eigen::VectorXd dual_solution =
			J_p(returnvalues.solution, PowerX);
	const double modulus_at_one = modul(1);
	const double _ANorm = ::pow(2, 1.+ 1./val_NormY); //_A.norm();
	while (!StopCriterion) {
		BOOST_LOG_TRIVIAL(debug)
				<< "#" << returnvalues.NumberOuterIterations
				<< " at " << returnvalues.solution.transpose()
				<< " in " << returnvalues.residual.transpose()
				<< " with residual of " << returnvalues.residuum;
		double alpha = 0.;
		if (dual_solution.isApproxToConstant(0, TolX)) {
			alpha =
					C * (::pow(val_DualNormX, val_NormX - 1.)
						/ ::pow(_ANorm,val_NormX))
					* ::pow(returnvalues.residuum, val_NormX - val_NormY);
		} else {
			const double xnorm = NormX(returnvalues.solution);
			alpha = C/(::pow(2., val_DualNormX) * G * _ANorm)
					* (returnvalues.residuum / xnorm);
			const double lambda = std::min(modulus_at_one, alpha);
			// find intermediate value in smoothness modulus to match lambda
			const double tau = calculateMatchingTau(modul, lambda);
			// calculate step width
			alpha = (tau/_ANorm)
					* (::pow( xnorm, val_NormX-1.)
						/ ::pow(returnvalues.residuum, val_NormY - 1.));
		}

		// iterate: J_p (x_{n+1})
		BOOST_LOG_TRIVIAL(trace)
			<< "Step width is " << alpha;
		BOOST_LOG_TRIVIAL(trace)
			<< "A^* is " << _A.transpose();
		BOOST_LOG_TRIVIAL(trace)
			<< "j_r (residual) is "
			<< j_r( returnvalues.residual, PowerY).transpose();
		const Eigen::VectorXd u =
				_A.transpose() * j_r( returnvalues.residual, PowerY);
		BOOST_LOG_TRIVIAL(trace)
			<< "u is " << u.transpose();
		dual_solution -= alpha * u;

		// finally map back from X^{\conj} to X: x_{n+1}
		returnvalues.solution =
				J_q(dual_solution, DualPowerX);

		// update residual
		returnvalues.residuum = calculateResidual(
				returnvalues.solution,
				_A,
				_y,
				NormY,
				returnvalues.residual);

		// check iterations count
		++returnvalues.NumberOuterIterations;
		StopCriterion =
				(returnvalues.NumberOuterIterations >= MaxOuterIterations)
				|| (fabs(returnvalues.residuum) <= TolY);
	}

	return returnvalues;
}
