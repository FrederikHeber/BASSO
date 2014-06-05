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
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <limits>

#include "BregmanFunctional.hpp"
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
		const Eigen::VectorXd &_solution,
		const unsigned int _outputsteps
		) :
	GeneralMinimizer(
			_NormX,
			_NormY,
			_PowerX,
			_PowerY,
			_Delta,
			_maxiter,
			_solution,
			_outputsteps
			),
	tau(1.1),
	N(2)
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
		const Eigen::VectorXd &_y
		) const
{
//	NoCols = _A.innerSize();
//	NoRows = _A.outerSize();

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

	BregmanDistance Delta_p(val_NormX);
	double old_distance = 0.;
	if (!solution.isZero()) {
		old_distance = Delta_p(returnvalues.solution, solution)
			+ 1e4*BASSOTOLERANCE; // make sure its larger
	}

	// start building up search space 'U' with the search vectors 'u'
	Eigen::MatrixXd U(_A.outerSize(),N);
	Eigen::VectorXd alphas(N);
	unsigned int index = 0;
	while (!StopCriterion) {
		BOOST_LOG_TRIVIAL(debug)
				<< "#" << returnvalues.NumberOuterIterations
				<< " with residual of " << returnvalues.residuum;
		// check that distance truely decreases
		if (!solution.isZero()) {
			const double new_distance = Delta_p(returnvalues.solution, solution);
			BOOST_LOG_TRIVIAL(debug)
				<< "Delta_p(x_" << returnvalues.NumberOuterIterations
				<< ",x) is "
				<< new_distance;
//			assert( old_distance > new_distance );
			old_distance = new_distance;
		}
		BOOST_LOG_TRIVIAL(trace)
				<< "x_n is " << returnvalues.solution.transpose();
		BOOST_LOG_TRIVIAL(trace)
				<< "R_n is " << returnvalues.residual.transpose();
		BOOST_LOG_TRIVIAL(trace)
			<< "j_r (residual) is "
			<< j_r( returnvalues.residual, PowerY).transpose();

		// u=A'*DualityMapping(w,NormY,PowerY,TolX);
		Eigen::VectorXd u = j_r(returnvalues.residual, PowerY);
		BOOST_LOG_TRIVIAL(trace)
			<< "u is " << u.transpose();

		// uNorm=norm(u,DualNormX);
		const double uNorm = DualNormX(u);
		BOOST_LOG_TRIVIAL(trace)
			<< "uNorm is " << uNorm;

		// alpha=u'*x
		const double alpha =
				u.transpose() * _y;
		BOOST_LOG_TRIVIAL(trace)
			<< "alpha is " << alpha;

		// add u to U and alpha to alphas
		U.col(index) = _A.transpose()*u;
		alphas(index) = alpha;
		index = (index + 1) % N;

		Eigen::VectorXd tmin(N);
		tmin.setZero();
		{
			// tmin=fminunc(@(t) BregmanFunctional(t,Jx,u,alpha+d,DualNormX,DualPowerX,TolX),t0,BregmanOptions);
			BregmanFunctional bregman(val_DualNormX, TolX);
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
			  iter++;
			  status = gsl_multimin_fdfminimizer_iterate (s);

			  if (status)
				break;

			  status = gsl_multimin_test_gradient (s->gradient, TolY);

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
		// x=DualityMapping(Jx-tmin*u,DualNormX,DualPowerX,TolX);
		dual_solution -= U*tmin;
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

	// and return solution
	return returnvalues;
}
