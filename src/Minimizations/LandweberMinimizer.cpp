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
		const double _C,
		const unsigned int _maxiter,
		const unsigned int _outputsteps
		) :
	val_NormX(_NormX),
	val_NormY(_NormY),
	val_DualNormX(val_NormX/(val_NormX - 1.)),
	PowerX(_PowerX),
	PowerY(_PowerY),
	DualPowerX(PowerX/(PowerX - 1.)),
	Delta(_Delta),
	MaxOuterIterations(_maxiter),
	TolX(1e-6),
	TolY(Delta),
	TolFun(1e-12),
	C(_C),
	outputsteps(_outputsteps),
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


/** Calculate residual \a _A * \a _x0 - \a _y in given norm \a _NormY.
 *
 * \param _x0 current iteration point
 * \param _A matrix of inverse problem
 * \param _y right-hand side
 * \param _residual residual vector, updated after call
 * \return norm of residual
 */
double LandweberMinimizer::calculateResidual(
		const Eigen::VectorXd &_x0,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y,
		Eigen::VectorXd &_residual
		) const
{
	_residual = _A * _x0 - _y;
	return NormY(_residual);
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
	const double _ANorm = ::pow(2, 1.+ 1./val_NormY); //_A.norm();
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

		// print each solution
		if ((outputsteps != 0) &&
				(returnvalues.NumberOuterIterations % outputsteps == 0)) {
			{
				std::stringstream solution_file;
				solution_file << "solution"
						<< (returnvalues.NumberOuterIterations / outputsteps) << ".m";
				using namespace MatrixIO;
				std::ofstream ost(solution_file.str().c_str());
				if (ost.good())
					ost << returnvalues.solution;
				else {
					std::cerr << "Failed to open " << solution_file.str() << std::endl;
				}
			}
			{
				std::stringstream solution_file;
				solution_file << "projected_solution"
						<< (returnvalues.NumberOuterIterations / outputsteps) << ".m";
				using namespace MatrixIO;
				std::ofstream ost(solution_file.str().c_str());
				if (ost.good())
					ost << _A * returnvalues.solution;
				else {
					std::cerr << "Failed to open " << solution_file.str() << std::endl;
				}
			}

		}
	}

	return returnvalues;
}
