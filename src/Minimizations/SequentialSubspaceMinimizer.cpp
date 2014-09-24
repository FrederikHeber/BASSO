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
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <limits>

#include "BregmanFunctional.hpp"
#include "Database/Database.hpp"
#include "Database/Table.hpp"
#include "DualityMapping.hpp"
#include "LpNorm.hpp"
#include "MinimizationExceptions.hpp"
#include "Minimizations/BregmanDistance.hpp"

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
	tau(1.1),
	N(2),
	MatrixVectorProduct_subspace(MatrixVectorProduct),
	ScalarVectorProduct_subspace(ScalarVectorProduct)
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

/** Structure containing all parameters to call BregmanFunctional functions.
 *
 * This is required to use function minimization that only allows
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
static double
func(const gsl_vector *x, void *adata)
{
	struct BregmanParameters *params =
			static_cast<BregmanParameters *>(adata);
	for (int i=0;i<params->t.innerSize();++i)
		params->t[i] = gsl_vector_get(x, i);
	const double returnvalue =
			(params->bregman)(
					params->t,
					params->x,
					params->U,
					params->alpha,
					params->q);
	BOOST_LOG_TRIVIAL(trace)
		<< "func() evaluates to " << returnvalue;
	return returnvalue;
}

/** Static function to wrap call to BregmanFunctional::operator()().
 *
 */
static void
jacf(const gsl_vector *x, void *adata, gsl_vector *g)
{
	struct BregmanParameters *params =
			static_cast<BregmanParameters *>(adata);
	// update where to evaluate
	for (int i=0;i<params->t.innerSize();++i)
		params->t[i] = gsl_vector_get(x, i);
	Eigen::VectorXd grad =
			(params->bregman).gradient(
					params->t,
					params->x,
					params->U,
					params->alpha,
					params->q);
	for (int i=0; i<grad.innerSize();++i)
		gsl_vector_set(g, i, grad[i]);
	BOOST_LOG_TRIVIAL(trace)
		<< "grad() evaluates to " << grad.transpose();
}

/** Static function to wrap call to BregmanFunctional::operator()().
 *
 */
static void
funcjacf(const gsl_vector *x, void *adata, double *f, gsl_vector *g)
{
	struct BregmanParameters *params =
			static_cast<BregmanParameters *>(adata);
	// update where to evaluate
	for (int i=0;i<params->t.innerSize();++i)
		params->t[i] = gsl_vector_get(x, i);
	*f =
			(params->bregman)(
					params->t,
					params->x,
					params->U,
					params->alpha,
					params->q);
	BOOST_LOG_TRIVIAL(trace)
		<< "func() evaluates to " << *f;
	Eigen::VectorXd grad =
			(params->bregman).gradient(
					params->t,
					params->x,
					params->U,
					params->alpha,
					params->q);
	for (int i=0; i<grad.innerSize();++i)
		gsl_vector_set(g, i, grad[i]);
	BOOST_LOG_TRIVIAL(trace)
		<< "grad() evaluates to " << grad.transpose();
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

	// start building up search space 'U' with the search vectors 'u'
	Eigen::MatrixXd U = Eigen::MatrixXd::Zero(_A.outerSize(),N);
	Eigen::VectorXd alphas = Eigen::VectorXd::Zero(N);
	unsigned int index = 0;
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
		U.col(index) = MatrixVectorProduct(A_transposed,Jw);
		//U.col(index) *= 1./NormX(U.col(index));
		alphas(index) = alpha;
		index = (index + 1) % N;

		Eigen::VectorXd tmin(N);
		tmin.setZero();
		{
			// tmin=fminunc(@(t) BregmanFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
			BregmanFunctional bregman(
					DualNormX,
					J_q,
					MatrixVectorProduct_subspace,
					ScalarVectorProduct_subspace);
			// we perform a function minimization
			// with respect to t starting at t0
			BregmanParameters params(
					bregman,
					tmin,
					dual_solution,
					U,
					alphas,
					DualPowerX);
			size_t iter = 0;
			int status;

			const gsl_multimin_fdfminimizer_type *T;
			gsl_multimin_fdfminimizer *s;

			gsl_vector *x;
			gsl_multimin_function_fdf my_func;

			my_func.n = N;
			my_func.f = &func;
			my_func.df = &jacf;
			my_func.fdf = &funcjacf;
			my_func.params = &params;

			/* Starting point, x = (0,0) */
			x = gsl_vector_alloc (N);
			for (unsigned int i=0;i<N;++i)
			  gsl_vector_set (x, i, tmin(i));

			T = gsl_multimin_fdfminimizer_vector_bfgs;
			s = gsl_multimin_fdfminimizer_alloc (T, N);

			gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, TolY);

			do
			{
			  ++iter;
			  status = gsl_multimin_fdfminimizer_iterate (s);

			  if (status)
				break;

			  status = gsl_multimin_test_gradient (s->gradient, TolFun);

//				  if (status == GSL_SUCCESS)
//					printf ("Minimum found at:\n");
//
//				  printf ("%5d %.5f %.5f %10.5f\n", iter,
//						  gsl_vector_get (s->x, 0),
//						  gsl_vector_get (s->x, 1),
//						  s->f);

			}
			while (status == GSL_CONTINUE && iter < 100);

			// place solution at tmin
			for (size_t i=0;i<N;++i)
			  tmin(i) = gsl_vector_get (s->x, i);

			gsl_multimin_fdfminimizer_free (s);
			gsl_vector_free (x);

			BOOST_LOG_TRIVIAL(trace)
				<< "tmin is " << tmin.transpose();
		}
		per_iteration_tuple.replace( "stepwidth", tmin.norm());
		// x=DualityMapping(Jx-tmin*u,DualNormX,DualPowerX,TolX);
		dual_solution -= MatrixVectorProduct(U,tmin);
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
