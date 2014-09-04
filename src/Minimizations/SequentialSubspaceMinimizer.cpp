/*
 * SequentialSubspaceMinimizer.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "SequentialSubspaceMinimizer.hpp"

#include <boost/chrono.hpp>
#include <boost/log/trivial.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <limits>

#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "DualityMapping.hpp"
#include "LpNorm.hpp"
#include "MinimizationExceptions.hpp"
#include "Minimizations/BregmanDistance.hpp"
#include "Minimizations/Functions/FunctionMinimizer.hpp"
#include "Minimizations/Functions/MinimizationFunctional.hpp"
#include "Minimizations/SequentialSubspaceMinimizer_HyperplaneProjection.hpp"

// instantiate required template functions
CONSTRUCT_FUNCTIONMINIMIZER(Eigen::VectorXd)

SequentialSubspaceMinimizer::SequentialSubspaceMinimizer(
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
	N(2),
	MatrixVectorProduct_subspace(MatrixVectorProduct),
	ScalarVectorProduct_subspace(ScalarVectorProduct)
{}

void SequentialSubspaceMinimizer::setN(
		const unsigned int _N
		)
{
	// check that number of search directions is greater than 1
	if (_N < 1)
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("N");
	const_cast<unsigned int&>(N) = _N;
}

SequentialSubspaceMinimizer::ReturnValues
SequentialSubspaceMinimizer::operator()(
		const Eigen::VectorXd &_x0,
		const Eigen::MatrixXd &_A,
		const Eigen::VectorXd &_y,
		const Eigen::VectorXd &_solution
		)
{
//	NoCols = _A.innerSize();
//	NoRows = _A.outerSize();
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

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

	/// -# calculate some values prior to loop
	// Jx=DualityMapping(x,NormX,PowerX,TolX);
	Eigen::VectorXd dual_solution =
			J_p(returnvalues.solution, PowerX);
	BOOST_LOG_TRIVIAL(trace)
		<< "Jx_0 is " << dual_solution.transpose();
//	const double modulus_at_one = modul(1);
//	const double _ANorm = ::pow(2, 1.+ 1./val_NormY); //_A.norm();
//	BOOST_LOG_TRIVIAL(trace)
//		<< "_ANorm " << _ANorm;

	if ((_solution.innerSize() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
			<< "true solution is " << _solution.transpose() << std::endl;
	} else {
		BOOST_LOG_TRIVIAL(info)
			<< "true solution is " << _solution.transpose() << std::endl;
	}

	BregmanDistance Delta_p(NormX, J_p, PowerX, ScalarVectorProduct);
	double old_distance = 0.;
	if (!_solution.isZero()) {
		old_distance = Delta_p(returnvalues.solution, _solution, PowerX)
			+ 1e4*BASSOTOLERANCE; // make sure its larger
		BOOST_LOG_TRIVIAL(debug)
				<< "Starting Bregman distance is " << old_distance;
	}

	// build data tuple for iteration information
	Table& per_iteration_table = database.addTable("per_iteration");
	Table::Tuple_t per_iteration_tuple;
	per_iteration_tuple.insert( std::make_pair("p", val_NormX));
	per_iteration_tuple.insert( std::make_pair("r", val_DualNormX));
	per_iteration_tuple.insert( std::make_pair("N", (int)N));
	per_iteration_tuple.insert( std::make_pair("dim", (int)_x0.innerSize()));
	per_iteration_tuple.insert( std::make_pair("iteration", (int)0));
	per_iteration_tuple.insert( std::make_pair("stepwidth", 0.));
	per_iteration_tuple.insert( std::make_pair("relative_residual", 0.));
	per_iteration_tuple.insert( std::make_pair("error", 0.));
	per_iteration_tuple.insert( std::make_pair("bregman_distance", 0.));

	// build data tuple for overall information
	Table& overall_table = database.addTable("overall");
	Table::Tuple_t overall_tuple;
	overall_tuple.insert( std::make_pair("p", val_NormX));
	overall_tuple.insert( std::make_pair("r", val_DualNormX));
	overall_tuple.insert( std::make_pair("N", (int)N));
	overall_tuple.insert( std::make_pair("dim", (int)_x0.innerSize()));
	overall_tuple.insert( std::make_pair("iterations", (int)0));
	overall_tuple.insert( std::make_pair("relative_residual", 0.));
	overall_tuple.insert( std::make_pair("runtime", 0.));
	overall_tuple.insert( std::make_pair("matrix_vector_products", (int)0));
	overall_tuple.insert( std::make_pair("vector_vector_products", (int)0));
	overall_tuple.insert( std::make_pair("matrix_vector_products_subspace", (int)0));
	overall_tuple.insert( std::make_pair("vector_vector_products_subspace", (int)0));

	/// -# check stopping criterion
	const Eigen::MatrixXd & A_transposed = _A.transpose();
	bool StopCriterion = false;
	StopCriterion = (fabs(returnvalues.residuum/_ynorm) <= TolY);

	// reset inner state of problem has changed
	if (state.getDimension() != _A.outerSize())
		state.set(_A.outerSize(), N);
	while (!StopCriterion) {
		per_iteration_tuple.replace( "iteration", (int)returnvalues.NumberOuterIterations);
		BOOST_LOG_TRIVIAL(debug)
				<< "#" << returnvalues.NumberOuterIterations
				<< " with residual of " << returnvalues.residuum;
		BOOST_LOG_TRIVIAL(debug)
			<< "#" << returnvalues.NumberOuterIterations << ": "
			<< "||Ax_n-y||/||y|| is " << returnvalues.residuum/_ynorm;
		per_iteration_tuple.replace( "relative_residual", returnvalues.residuum/_ynorm);
		// check that distance truly decreases
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

		// Jw=DualityMapping(w,NormY,PowerY,TolX);
		Eigen::VectorXd Jw = j_r(returnvalues.residual, PowerY);
		BOOST_LOG_TRIVIAL(trace)
			<< "Jw= j_r (R_n) is " << Jw.transpose();

		// JwNorm=norm(w,DualNormX);
//		const double JwNorm = DualNormX(Jw);
//		BOOST_LOG_TRIVIAL(trace)
//			<< "wNorm is " << wNorm;

		// alpha=Jw'*y
		const double alpha =
				Jw.transpose() * _y;
		BOOST_LOG_TRIVIAL(trace)
			<< "alpha is " << alpha;

		// add u to U and alpha to alphas
		state.U.col(state.index) = MatrixVectorProduct(A_transposed,Jw);
		//U.col(index) *= 1./NormX(U.col(index));
		state.alphas(state.index) = alpha;
		state.index = (state.index + 1) % N;

		Eigen::VectorXd tmin(N);
		tmin.setZero();
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

			FunctionMinimizer<Eigen::VectorXd> minimizer(
				functional, tmin);

			tmin = minimizer(N, TolFun, tmin);

			BOOST_LOG_TRIVIAL(trace)
				<< "tmin is " << tmin.transpose();
		}
		per_iteration_tuple.replace( "stepwidth", tmin.norm());
		// x=DualityMapping(Jx-tmin*u,DualNormX,DualPowerX,TolX);
		dual_solution -= MatrixVectorProduct(state.U,tmin);
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
	overall_tuple.replace( "vector_vector_products", (int)ScalarVectorProduct.getCount() );
	overall_tuple.replace( "matrix_vector_products_subspace", (int)MatrixVectorProduct_subspace.getCount() );
	overall_tuple.replace( "vector_vector_products_subspace", (int)ScalarVectorProduct_subspace.getCount() );
	overall_table.addTuple(overall_tuple);

	// and return solution
	return returnvalues;
}

void SequentialSubspaceMinimizer::IterationState::set(
		const unsigned int _dimension,
		const unsigned int _N
		)
{
	U = Eigen::MatrixXd::Zero(_dimension,_N);
	alphas = Eigen::VectorXd::Zero(_N);
	index = 0;
	isInitialized = true;
}
