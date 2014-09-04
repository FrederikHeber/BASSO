/*
 * FunctionMinimizer_impl.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef FUNCTIONMINIMIZER_IMPL_HPP_
#define FUNCTIONMINIMIZER_IMPL_HPP_

#include "BassoConfig.h"

#include "Minimizations/Functions/FunctionMinimizer.hpp"

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

template <class T>
const T FunctionMinimizer<T>::operator()(
		const unsigned int _N,
		const double _Tol,
		const T &_startvalue)
{
	size_t iter = 0;
	int status;

	const gsl_multimin_fdfminimizer_type *minimizer;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *x;
	gsl_multimin_function_fdf my_func;

	my_func.n = _N;
	my_func.f = &FunctionMinimizer_FunctionCaller<T>;
	my_func.df = &FunctionMinimizer_GradientCaller<T>;
	my_func.fdf = &FunctionMinimizer_FunctionGradientCaller<T>;
	my_func.params = this;

	/* Starting point, x = (0,0) */
	x = gsl_vector_alloc (_N);
	functional.convertFromInternalType(_startvalue, x);

	minimizer = gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc (minimizer, _N);

	gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, _Tol);

	do
	{
	  ++iter;
	  status = gsl_multimin_fdfminimizer_iterate (s);

	  if (status)
		break;

	  status = gsl_multimin_test_gradient (s->gradient, _Tol);

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
	functional.convertToInternalType(value, s->x);

	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);

	return value;
}


#endif /* FUNCTIONMINIMIZER_IMPL_HPP_ */
