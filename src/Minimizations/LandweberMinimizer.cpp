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
#include <fstream>
#include <sstream>

#include "MatrixIO/MatrixIO.hpp"
#include "MinimizationExceptions.hpp"

LandweberMinimizer::LandweberMinimizer(
		const double _NormX,
		const double _NormY,
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
	C(0.9),
	modul(val_DualNormX)
{}

void LandweberMinimizer::setC(const double _C)
{
	if ((_C <= 0.)) // || ( _C > 1.))
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("C");

	const_cast<double&>(C) = _C;
}

GeneralMinimizer::ReturnValues
LandweberMinimizer::operator()(
		const Eigen::VectorXd &_x0,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y
		) const
{
	if ((_A.innerSize() < 10) && (_A.outerSize() < 10)) {
		BOOST_LOG_TRIVIAL(info) << "Calculating "
				<< _A << "*" << _x0.transpose() << "-" << _y.transpose();
	} else {
		BOOST_LOG_TRIVIAL(trace) << "Calculating "
				<< _A << "*" << _x0.transpose() << "-" << _y.transpose();
	}

	// G constant used in theoretical step width
	double G;
	if (val_DualNormX < 2.) {
		G = ::pow(2., 2. - val_DualNormX);
	} else {
//		PowerX = 2.;
//		DualPowerX = 2.;
		G = val_DualNormX - 1.;
	}
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
	Eigen::VectorXd dual_solution =
			J_p(returnvalues.solution, PowerX);
	const double modulus_at_one = modul(1);
	const double _ANorm = _A.norm(); //::pow(2, 1.+ 1./val_NormY);
	BOOST_LOG_TRIVIAL(trace)
		<< "_ANorm " << _ANorm;

	/// -# loop over stopping criterion
	while (!StopCriterion) {
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
		const Eigen::VectorXd u =
				_A.transpose() * j_r( returnvalues.residual, PowerY);
		BOOST_LOG_TRIVIAL(trace)
			<< "u is " << u.transpose();
		double alpha = 0.;
		// use step width used in theoretical proof
		// (F. Schöpfer, 11.4.2014) too conservative! Line search instead
		if (0) {
		if (dual_solution.isApproxToConstant(0, TolX)) {
			const double q_p = ::pow(val_DualNormX, val_NormX - 1.);
			BOOST_LOG_TRIVIAL(trace)
				<< "q_p " << q_p;
			const double A_p = ::pow(_ANorm,val_NormX);
			BOOST_LOG_TRIVIAL(trace)
				<< "A_p " << A_p;
			double R_p = 0.;
			if (val_NormX == LpNorm::Infinity) {
				R_p = returnvalues.residual.array().abs().maxCoeff();
			} else {
				R_p = ::pow(returnvalues.residuum, val_NormX);
			}
			double R_r = 0.;
			if (val_NormY == LpNorm::Infinity) {
				R_r = returnvalues.residual.array().abs().maxCoeff();
			} else {
				R_r = ::pow(returnvalues.residuum, val_NormY);
			}
			BOOST_LOG_TRIVIAL(trace)
				<< "R_p " << R_p << ", R_r " << R_r;
			alpha = // C
					0.9 * (q_p / A_p) * (R_p/R_r);
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
			double R_r = 0.;
			if (val_NormY == LpNorm::Infinity) {
				R_r = returnvalues.residual.array().abs().maxCoeff()
						/ returnvalues.residuum;
			} else {
				R_r = ::pow(returnvalues.residuum, val_NormY - 1.);
			}
			BOOST_LOG_TRIVIAL(trace)
				<< "R_r " << R_r;
			alpha = (tau/_ANorm) * (x_p / R_r);
		}
		}
		// get alpha from line search
		alpha = calculateOptimalStepwidth(
				 returnvalues.solution,
				 u,
				 _A,
				 _y,
				 alpha);

		// iterate: J_p (x_{n+1})
		BOOST_LOG_TRIVIAL(trace)
			<< "Step width is " << alpha;
		dual_solution -= alpha * u;
		BOOST_LOG_TRIVIAL(trace)
				<< "x^*_n+1 is " << dual_solution.transpose();

		// finally map back from X^{\conj} to X: x_{n+1}
		returnvalues.solution =
				J_q(dual_solution, DualPowerX);
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

	return returnvalues;
}
