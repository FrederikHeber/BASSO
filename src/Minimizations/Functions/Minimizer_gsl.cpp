/*
 * Minimizer_gsl.cpp
 *
 *  Created on: Nov 4, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Minimizer_gsl.hpp"

#include <sstream>

#include "Log/Logging.hpp"

Minimizer<gsl_vector>::Minimizer(
		const unsigned int _N
		) :
		N(_N),
		maxiterations(100),
		optimum(0),
		tempoptimum(N),
		tempgradient(N)
{
	my_func.n = N;
	my_func.f = &Minimizer<gsl_vector>::FunctionCaller;
	my_func.df = &Minimizer<gsl_vector>::GradientCaller;
	my_func.fdf = &Minimizer<gsl_vector>::FunctionGradientCaller;
	my_func.params = this;

	minimizer = gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc (minimizer, N);
}

Minimizer<gsl_vector>::~Minimizer()
{
	gsl_multimin_fdfminimizer_free (s);
}

enum Minimization::GradientStatus
Minimizer<gsl_vector>::checkGradient(
		const double _Tol) const
{
	int status = gsl_multimin_test_gradient (
			gsl_multimin_fdfminimizer_gradient(s),
			_Tol);
	switch (status) {
		case GSL_SUCCESS:
			return Minimization::gradient_success;
			break;
		default:
		case GSL_CONTINUE:
			return Minimization::gradient_continue;
			break;
		case GSL_ENOPROG:
			return Minimization::error_noprogress;
			break;
	}
}

const unsigned int
Minimizer<gsl_vector>::minimize(
		const double _Tol,
		array_type &_startvalue,
		const check_function_t &_checkfunction
		)
{
	unsigned int iter = 0;
	enum Minimization::GradientStatus status;

	/* Starting point, x = (0,0) */
	gsl_vector *x;
	x = gsl_vector_alloc (N);
	convertArrayTypeToInternalType(_startvalue, x);

	gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, _Tol);

	do
	{
		++iter;
		int gsl_status = gsl_multimin_fdfminimizer_iterate (s);

		optimum = gsl_multimin_fdfminimizer_minimum(s);
		convertInternalTypeToArrayType(
				gsl_multimin_fdfminimizer_x(s), tempoptimum);
		convertInternalTypeToArrayType(
				gsl_multimin_fdfminimizer_gradient(s), tempgradient);

		{
			std::stringstream iterate;
			const gsl_vector * const currentiterate =
					gsl_multimin_fdfminimizer_x(s);
			iterate << "Current iterate #" << iter << ":";
			for (unsigned int i=0;i<N;++i)
				iterate << " " << gsl_vector_get(currentiterate, i);
			BOOST_LOG_TRIVIAL(debug)
				<< iterate.str();
		}

		if (gsl_status == GSL_ENOPROG) {
			BOOST_LOG_TRIVIAL(warning)
					<< "gsl_multimin could not improve iterate anymore.";
			break;
		}

		status = _checkfunction(_Tol);

	}
	while ((status == Minimization::gradient_continue) && iter < maxiterations);

	BOOST_LOG_TRIVIAL(debug)
		<< "Inner iteration took " << iter << " steps";

	// place solution at tmin
	convertInternalTypeToArrayType(
			gsl_multimin_fdfminimizer_x(s), _startvalue);

	gsl_vector_free (x);

	return iter;
}

double
Minimizer<gsl_vector>::FunctionCaller(
		const gsl_vector *x,
		void *adata)
{
	// obtain minimizer object
	struct Minimizer<gsl_vector> *minimizer =
			static_cast<Minimizer<gsl_vector> *>(adata);
	// convert given value to something interpretable
	minimizer->convertInternalTypeToArrayType(x, minimizer->tempoptimum);
	// evaluate the function (with array_type)
	return minimizer->function_evaluator(minimizer->tempoptimum);
}

void
Minimizer<gsl_vector>::GradientCaller(
		const gsl_vector *x,
		void *adata,
		gsl_vector *g)
{
	// obtain minimizer object
	struct Minimizer<gsl_vector> *minimizer =
			static_cast<Minimizer<gsl_vector> *>(adata);
	// convert given value to something interpretable
	minimizer->convertInternalTypeToArrayType(x, minimizer->tempoptimum);
	// evaluate the gradient
	minimizer->tempgradient =
			minimizer->gradient_evaluator(minimizer->tempoptimum);
	// convert gradient and store in g
	minimizer->convertArrayTypeToInternalType(minimizer->tempgradient, g);
}

void
Minimizer<gsl_vector>::FunctionGradientCaller(
		const gsl_vector *x,
		void *adata,
		double *f,
		gsl_vector *g)
{
	// obtain minimizer object
	struct Minimizer<gsl_vector> *minimizer =
			static_cast<Minimizer<gsl_vector> *>(adata);
	// convert given value to something interpretable
	minimizer->convertInternalTypeToArrayType(x, minimizer->tempoptimum);
	// evaluate the function and store in f
	*f = minimizer->function_evaluator(minimizer->tempoptimum);
	// evaluate the gradient
	minimizer->tempgradient =
			minimizer->gradient_evaluator(minimizer->tempoptimum);
	// convert gradient and store in g
	minimizer->convertArrayTypeToInternalType(minimizer->tempgradient, g);
}

void
Minimizer<gsl_vector>::convertArrayTypeToInternalType(
		const array_type & _t,
		gsl_vector * const _x)
{
	for(unsigned int i=0;i<_x->size;++i)
		gsl_vector_set(_x,i, _t[i]);
}

void
Minimizer<gsl_vector>::convertInternalTypeToArrayType(
		const gsl_vector * const _x,
		array_type &_t)
{
	for(unsigned int i=0;i<_x->size;++i)
		_t[i] = gsl_vector_get(_x,i);
}
