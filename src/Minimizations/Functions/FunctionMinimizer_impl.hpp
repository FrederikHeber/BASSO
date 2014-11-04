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

#include <boost/log/trivial.hpp>

static
int checkGradient(
		const gsl_multimin_fdfminimizer * const _s,
		const double _Tol)
{
	return gsl_multimin_test_gradient (_s->gradient, _Tol);
}

template <class T>
int FunctionMinimizer<T>::checkWolfeConditions(
		const double _startvalue,
		const gsl_vector * const _startgradient,
		const Wolfe_indexset_t &_Wolfe_indexset,
		const gsl_multimin_fdfminimizer * const _s) const
{
	bool conditions_fulfilled = true;
	assert( !_Wolfe_indexset.empty() );
	for (Wolfe_indexset_t::const_iterator iter = _Wolfe_indexset.begin();
			iter != _Wolfe_indexset.end(); ++iter) {
		const double component = gsl_vector_get(_s->x, *iter);
		// 1. positivity of step width component
		BOOST_LOG_TRIVIAL(trace)
			<< "1. Positivity: " << component << " > " << constant_positivity;
		conditions_fulfilled &=
				component > constant_positivity;
		// 2. still descent in the current gradient component
		BOOST_LOG_TRIVIAL(trace)
			<< "2. Still descent: " << gsl_vector_get(_s->gradient, *iter)
			<< " < " << std::numeric_limits<double>::epsilon();
		conditions_fulfilled &=
				(gsl_vector_get(_s->gradient, *iter)
						< std::numeric_limits<double>::epsilon());
		// 3. stronger descent than linear interpolation
		const double linearinterpolate = _startvalue
				+ constant_interpolation * component *
				gsl_vector_get(_startgradient, *iter);
		BOOST_LOG_TRIVIAL(trace)
			<< "3. Stronger than linear: " << _s->f
			<< " < " << linearinterpolate;
		conditions_fulfilled &=
				(_s->f < linearinterpolate);
	}
	if (conditions_fulfilled)
		return GSL_SUCCESS;
	else
		return GSL_CONTINUE;
}

template <class T>
const unsigned int
FunctionMinimizer<T>::operator()(
		const unsigned int _N,
		const double _Tol,
		T &_startvalue)
{
	check_function_t checkfunction =
			boost::bind(&checkGradient,
					_1, _2);
	return performMinimization(
			_N,_Tol,_startvalue,
			checkfunction);
}

template <class T>
const unsigned int
FunctionMinimizer<T>::operator()(
		const unsigned int _N,
		const double _Tol,
		const Wolfe_indexset_t &_Wolfe_indexset,
		T &_startvalue)
{
	gsl_vector *x = gsl_vector_alloc (_N);
	const double functionzero = FunctionMinimizer_FunctionCaller<T>(x, this);
	gsl_vector *gradientzero = gsl_vector_alloc (_N);
	FunctionMinimizer_GradientCaller<T>(x, this, gradientzero);
	gsl_vector_free(x);

	check_function_t checkfunction =
			boost::bind(&FunctionMinimizer<T>::checkWolfeConditions,
					boost::cref(*this),
					functionzero,
					gradientzero,
					boost::cref(_Wolfe_indexset),
					_1);
	const unsigned int iterations = performMinimization(
			_N,_Tol,_startvalue,
			checkfunction);

	gsl_vector_free(gradientzero);

	return iterations;
}

template <class T>
const unsigned int
FunctionMinimizer<T>::performMinimization(
		const unsigned int _N,
		const double _Tol,
		T &_startvalue,
		const check_function_t &_checkfunction
		)
{
	unsigned int iter = 0;
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
	std::vector<double> tempvector(_N);
	functional.convertToArrayType(_startvalue, tempvector);
	convertToInternalType(tempvector, x);

	minimizer = gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc (minimizer, _N);

	gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, _Tol);

	do
	{
		++iter;
		status = gsl_multimin_fdfminimizer_iterate (s);

		{
			std::stringstream iterate;
			const gsl_vector * const currentiterate =
					gsl_multimin_fdfminimizer_x(s);
			iterate << "Current iterate #" << iter << ":";
			for (unsigned int i=0;i<_N;++i)
				iterate << " " << gsl_vector_get(currentiterate, i);
			BOOST_LOG_TRIVIAL(debug)
				<< iterate.str();
		}

		if (status == GSL_ENOPROG) {
			BOOST_LOG_TRIVIAL(warning)
					<< "gsl_multimin could not improve iterate anymore.";
			break;
		}

		status = _checkfunction(s, _Tol);

//				  if (status == GSL_SUCCESS)
//					printf ("Minimum found at:\n");
//
//				  printf ("%5d %.5f %.5f %10.5f\n", iter,
//						  gsl_vector_get (s->x, 0),
//						  gsl_vector_get (s->x, 1),
//						  s->f);

	}
	while (status == GSL_CONTINUE && iter < maxiterations);

	BOOST_LOG_TRIVIAL(debug)
		<< "Inner iteration took " << iter << " steps";

	// place solution at tmin
	convertFromInternalType(s->x, tempvector);
	functional.convertFromArrayType(tempvector, _startvalue);

	gsl_multimin_fdfminimizer_free (s);
	gsl_vector_free (x);

	return iter;
}

template <class T>
void
FunctionMinimizer<T>::convertToInternalType(
		const array_type & _t,
		gsl_vector * const _x) const
{
	for(unsigned int i=0;i<_x->size;++i)
		gsl_vector_set(_x,i, _t[i]);
}

template <class T>
void
FunctionMinimizer<T>::convertFromInternalType(
		const gsl_vector * const _x,
		array_type &_t) const
{
	for(unsigned int i=0;i<_x->size;++i)
		_t[i] = gsl_vector_get(_x,i);
}


#endif /* FUNCTIONMINIMIZER_IMPL_HPP_ */
