/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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
#include "Minimizations/Functions/Minimizers/MinimizerExceptions.hpp"


template<class T>
inline boost::log::formatting_ostream&
operator<<(boost::log::formatting_ostream&ost, const gsl_vector * values)
{
	for (size_t i=0;i<values->size;++i)
		ost << gsl_vector_get(values, i) << " ";
	return ost;
}

Minimizer<gsl_vector>::Minimizer(
		const unsigned int _N
		) :
		FunctionMinimizer(_N),
		N(_N)
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

enum FunctionMinimizer::GradientStatus
Minimizer<gsl_vector>::checkGradient(
		const double _Tol) const
{
	int status = gsl_multimin_test_gradient (
			gsl_multimin_fdfminimizer_gradient(s),
			_Tol);
	switch (status) {
		case GSL_SUCCESS:
			return FunctionMinimizer::gradient_success;
			break;
		default:
		case GSL_CONTINUE:
			return FunctionMinimizer::gradient_continue;
			break;
		case GSL_ENOPROG:
			return FunctionMinimizer::error_noprogress;
			break;
	}
}

static void doWarnIterate()
{
	static bool repeating = false;
	if (!repeating) {
		LOG(warning, "gsl_multimin could not improve iterate anymore, not warning any longer.");
		repeating = true;
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
	enum FunctionMinimizer::GradientStatus status;

	/* Starting point, x = (0,0) */
	gsl_vector *x;
	x = gsl_vector_alloc (N);
	convertArrayTypeToInternalType(_startvalue, x);

	// use inexact line searches here with 0.1 precision
	gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 0.1);

	do
	{
		++iter;
		int gsl_status = gsl_multimin_fdfminimizer_iterate (s);

		optimum = gsl_multimin_fdfminimizer_minimum(s);
		if (isnan(optimum) || isinf(optimum))
			throw MinimizerIllegalNumber_exception()
			<< MinimizerIllegalNumber_variablename("optimum");
		convertInternalTypeToArrayType(
				gsl_multimin_fdfminimizer_x(s), tempoptimum);
		convertInternalTypeToArrayType(
				gsl_multimin_fdfminimizer_gradient(s), tempgradient);

		{
			const gsl_vector * const currentiterate =
					gsl_multimin_fdfminimizer_x(s);
			for (unsigned int i=0;i<N;++i) {
				const double &value = gsl_vector_get(currentiterate, i);
				if (isnan(value) || isinf(value)) {
					std::stringstream iterate;
					iterate << "currentiterate, #" << i;
					throw MinimizerIllegalNumber_exception()
					<< MinimizerIllegalNumber_variablename(iterate.str());
				}
			}
		}
		LOG(trace, "Current iterate #" << iter << ":" << " " << gsl_multimin_fdfminimizer_x(s));

		if (gsl_status == GSL_ENOPROG) {
			doWarnIterate();
			break;
		}

		status = _checkfunction(_Tol);

	}
	while ((status == FunctionMinimizer::gradient_continue) && iter < maxiterations);

	LOG(trace, "Inner iteration took " << iter << " steps");

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
