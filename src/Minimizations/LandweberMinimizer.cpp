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

#include "MinimizationExceptions.hpp"

LandweberMinimizer::LandweberMinimizer(
		const double _NormX,
		const double _NormY,
		const double _PowerY,
		const double _Delta,
		const double _C,
		const unsigned int _maxiter
		) :
	val_NormX(_NormX),
	val_NormY(_NormY),
	val_DualNormX(val_NormX/(val_NormX - 1.)),
	PowerX(val_NormX),
	PowerY(_PowerY),
	DualPowerX(val_DualNormX),  // PowerX/(PowerX - 1.)
	Delta(_Delta),
	MaxOuterIterations(_maxiter),
	TolX(1e-6),
	TolY(Delta),
	TolFun(1e-12),
	C(_C),
	NormX(val_NormX),
	NormY(val_NormY),
	DualNormX(val_DualNormX),
	J_p(val_NormX),
	J_q(val_DualNormX),
	j_r(val_NormY),
	modul(val_DualNormX)
{
	if ((C <= 0.)) // || ( C > 1.))
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("C");

	BOOST_LOG_TRIVIAL(debug)
		<< "p is " << val_NormX
		<< ", q is " << val_DualNormX
		<< ", r is " << val_NormY
		<< ", power of J_p is " <<  PowerX
		<< ", power of J_q is " <<  DualPowerX
		<< ", power of J_r is " <<  PowerY;

	// set tolerances values
	J_p.setTolerance(TolX);
	J_q.setTolerance(TolX);
	j_r.setTolerance(TolY);
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
	BOOST_LOG_TRIVIAL(debug) << "Calculating "
			<< _A << "*" << _x0 << "-" << _y;
	BOOST_LOG_TRIVIAL(trace)
		<< "A^* is " << _A.transpose();

	// G constant used in theoretical step width
	double G;
	if (val_DualNormX < 2.) {
		G = ::pow(2., 2. - val_DualNormX);
	} else {
//		PowerX = 2.;
//		DualPowerX = 2.;
		G = val_DualNormX - 1.;
	}

	/// -# initialize return structure
	ReturnValues returnvalues;
	returnvalues.NumberOuterIterations = 0;
	// set iterate 'x' as start vector 'x0'
	returnvalues.solution = _x0;
	// calculate starting residual and norm
	returnvalues.residuum = calculateResidual(
			_x0, _A, _y, NormY,
			returnvalues.residual);

	/// -# check stopping criterion
	bool StopCriterion = false;
	StopCriterion = (fabs(returnvalues.residuum) <= TolY);

	// calculate some values prior to loop
	Eigen::VectorXd dual_solution =
			J_p(returnvalues.solution, PowerX);
	const double modulus_at_one = modul(1);
	const double _ANorm = ::pow(2, 1.+ 1./val_NormY); //_A.norm();
	BOOST_LOG_TRIVIAL(trace)
		<< "_ANorm " << _ANorm;

	/// -# loop over stopping criterion
	while (!StopCriterion) {
		BOOST_LOG_TRIVIAL(debug)
				<< "#" << returnvalues.NumberOuterIterations
				<< " at " << returnvalues.solution.transpose()
				<< " in " << returnvalues.residual.transpose()
				<< " with residual of " << returnvalues.residuum;
		double alpha = 0.;
		if (dual_solution.isApproxToConstant(0, TolX)) {
			const double q_p = ::pow(val_DualNormX, val_NormX - 1.);
			BOOST_LOG_TRIVIAL(trace)
				<< "q_p " << q_p;
			const double A_p = ::pow(_ANorm,val_NormX);
			BOOST_LOG_TRIVIAL(trace)
				<< "A_p " << A_p;
			const double R_p_r = ::pow(returnvalues.residuum, val_NormX - val_NormY);
			BOOST_LOG_TRIVIAL(trace)
				<< "R_p_r " << R_p_r;
			alpha = // C
					0.9 * (q_p / A_p) * R_p_r;
		} else {
			const double xnorm = NormX(returnvalues.solution);
			BOOST_LOG_TRIVIAL(trace)
				<< "xnorm " << xnorm;
			const double two_q = ::pow(2., val_DualNormX);
			BOOST_LOG_TRIVIAL(trace)
				<< "two_q " << two_q;
			alpha = C/(two_q * G * _ANorm)
					* (returnvalues.residuum / xnorm);
			BOOST_LOG_TRIVIAL(trace)
				<< "initial lambda " << alpha;
			const double lambda = std::min(modulus_at_one, alpha);
			// find intermediate value in smoothness modulus to match lambda
			const double tau = calculateMatchingTau(modul, lambda);
			// calculate step width
			const double x_p = ::pow( xnorm, val_NormX-1.);
			BOOST_LOG_TRIVIAL(trace)
				<< "x_p " << x_p;
			const double R_r = ::pow(returnvalues.residuum, val_NormY - 1.);
			BOOST_LOG_TRIVIAL(trace)
				<< "R_r " << R_r;
			alpha = (tau/_ANorm) * (x_p / R_r);
		}

		// iterate: J_p (x_{n+1})
		BOOST_LOG_TRIVIAL(trace)
			<< "Step width is " << alpha;
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
