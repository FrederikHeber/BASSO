/*
 * LandweberMinimizer.cpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#include "LandweberMinimizer.hpp"

#include <algorithm>
#include <boost/chrono.hpp>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <limits>
#include <fstream>
#include <sstream>

#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "MinimizationExceptions.hpp"
#include "Minimizations/BregmanDistance.hpp"

LandweberMinimizer::LandweberMinimizer(
		const double _NormX,
		const double _NormY,
		const double _PowerX,
		const double _PowerY,
		const double _Delta,
		const unsigned int _maxiter,
		Database &_database,
		const unsigned int _outputsteps
		) :
	GeneralMinimizer(
				_NormX,
				_NormY,
				_PowerX,
				_PowerY,
				_Delta,
				_maxiter,
				_database,
				_outputsteps
				),
	C(0.9),
	modul(val_DualNormX),
	useOptimalStepwidth(true)
{}

void LandweberMinimizer::setuseOptimalStepwidth(
	const bool _useOptimalStepwidth)
{
	const_cast<bool &>(useOptimalStepwidth) = _useOptimalStepwidth;
}

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
		const Eigen::VectorXd &_y,
		const Eigen::VectorXd &_solution
		)
{
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

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
	const double _ynorm = NormY(_y);

	// build data tuple for iteration information
	Table& per_iteration_table = database.addTable("per_iteration");
	Table::Tuple_t per_iteration_tuple;
	per_iteration_tuple.insert( std::make_pair("p", val_NormX));
	per_iteration_tuple.insert( std::make_pair("r", val_DualNormX));
	per_iteration_tuple.insert( std::make_pair("dim", (int)_x0.innerSize()));
	per_iteration_tuple.insert( std::make_pair("iteration", (int)0));
	per_iteration_tuple.insert( std::make_pair("stepwidth", (int)0));
	per_iteration_tuple.insert( std::make_pair("relative_residual", 0.));
	per_iteration_tuple.insert( std::make_pair("error", 0.));
	per_iteration_tuple.insert( std::make_pair("bregman_distance", 0.));

	// build data tuple for overall information
	Table& overall_table = database.addTable("overall");
	Table::Tuple_t overall_tuple;
	overall_tuple.insert( std::make_pair("p", val_NormX));
	overall_tuple.insert( std::make_pair("r", val_DualNormX));
	overall_tuple.insert( std::make_pair("dim", (int)_x0.innerSize()));
	overall_tuple.insert( std::make_pair("iterations", (int)0));
	overall_tuple.insert( std::make_pair("relative_residual", 0.));
	overall_tuple.insert( std::make_pair("runtime", 0.));
	overall_tuple.insert( std::make_pair("matrix_vector_products", (int)0));
	overall_tuple.insert( std::make_pair("runtime_matrix_vector_products", 0.));
	overall_tuple.insert( std::make_pair("vector_vector_products", (int)0));
	overall_tuple.insert( std::make_pair("runtime_vector_vector_products", 0.));

	/// -# check stopping criterion
	bool StopCriterion = false;
	StopCriterion = (fabs(returnvalues.residuum/_ynorm) <= TolY);

	// calculate some values prior to loop
	Eigen::VectorXd dual_solution =
			J_p(returnvalues.solution, PowerX);
	const double modulus_at_one = modul(1);
	double _ANorm = 0.;
	if (!useOptimalStepwidth) {
		_ANorm = _A.norm(); //::pow(2, 1.+ 1./val_NormY);
		BOOST_LOG_TRIVIAL(trace)
			<< "_ANorm " << _ANorm;
	}
	BregmanDistance Delta_p(NormX, J_p, val_NormX, ScalarVectorProduct);
	double old_distance = 0.;
	if (!_solution.isZero()) {
		old_distance = Delta_p(returnvalues.solution, _solution, PowerX)
			+ 1e4*BASSOTOLERANCE; // make sure its larger
		BOOST_LOG_TRIVIAL(debug)
				<< "Starting Bregman distance is " << old_distance;
	}

	/// -# loop over stopping criterion
	while (!StopCriterion) {
		per_iteration_tuple.replace( "iteration", (int)returnvalues.NumberOuterIterations);
		BOOST_LOG_TRIVIAL(debug)
				<< "#" << returnvalues.NumberOuterIterations
				<< " with residual of " << returnvalues.residuum;
		BOOST_LOG_TRIVIAL(debug)
			<< "#" << returnvalues.NumberOuterIterations << ": "
			<< "||Ax_n-y||/||y|| is " << returnvalues.residuum/_ynorm;
		per_iteration_tuple.replace( "relative_residual", returnvalues.residuum/_ynorm);
		// check that distance truely decreases
		if (!_solution.isZero()) {
			const double new_distance =
					Delta_p(returnvalues.solution, _solution, PowerX);
			BOOST_LOG_TRIVIAL(debug)
				<< "#" << returnvalues.NumberOuterIterations << ": "
				<< "Delta_p(x_n,x) is "
				<< new_distance;
			per_iteration_tuple.replace( "bregman_distance", new_distance);
//			assert( old_distance > new_distance );
			old_distance = new_distance;
			BOOST_LOG_TRIVIAL(debug)
				<< "#" << returnvalues.NumberOuterIterations << ": "
				<< "||x_n-x|| is " << NormX(returnvalues.solution-_solution);
			per_iteration_tuple.replace( "error", NormX(returnvalues.solution-_solution));
		}
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
		if (!useOptimalStepwidth) {
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
		} else {
			// get alpha from line search
			alpha = calculateOptimalStepwidth(
					 returnvalues.solution,
					 u,
					 _A,
					 _y,
					 alpha);
		}

		per_iteration_tuple.replace( "stepwidth", alpha);

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
				|| (fabs(returnvalues.residuum/_ynorm) <= TolY);

		// print intermediat solution
		printIntermediateSolution(
				returnvalues.solution,
				_A,
				returnvalues.NumberOuterIterations);

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
	overall_tuple.replace( "runtime_matrix_vector_products", MatrixVectorProduct.getTiming() );
	overall_tuple.replace( "vector_vector_products", (int)ScalarVectorProduct.getCount() );
	overall_tuple.replace( "runtime_vector_vector_products", ScalarVectorProduct.getTiming() );
	overall_table.addTuple(overall_tuple);

	return returnvalues;
}
