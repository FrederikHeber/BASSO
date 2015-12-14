/*
 * FunctionalMinimizer_exactLinesearch_impl.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef FUNCTIONALMINIMIZER_EXACTLINESEARCH_IMPL_HPP_
#define FUNCTIONALMINIMIZER_EXACTLINESEARCH_IMPL_HPP_

#include "BassoConfig.h"

#include "Math/Helpers.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer_exactLinesearch.hpp"

#include <boost/bind.hpp>

#include <boost/log/trivial.hpp>

template <class S>
FunctionalMinimizer_exactLinesearch<S>::FunctionalMinimizer_exactLinesearch(
		const MinimizationFunctional<S> &_functional,
		FunctionMinimizer::ptr_t &_minimizer,
		S &_value) :
		FunctionalMinimizer<S>::FunctionalMinimizer(
				_functional, _minimizer, _value)
{}

template <class S>
const unsigned int
FunctionalMinimizer_exactLinesearch<S>::operator()(
		const unsigned int _N,
		const double _Tol,
		S &_startvalue) const
{
	// use the default check to stop line search
	check_function_t checkfunction =
			boost::bind(&FunctionMinimizer::checkGradient,
					boost::cref(FunctionalMinimizer<S>::minimizer), _1);

	array_type functionargument(_N, 0.);
	functional.convertInternalTypeToArrayType(_startvalue, functionargument);
	const unsigned int iterations =
			minimizer->minimize(_Tol, functionargument, checkfunction);
	functional.convertArrayTypeToInternalType(functionargument, _startvalue);
	return iterations;
}

#endif /* FUNCTIONALMINIMIZER_EXACTLINESEARCH_IMPL_HPP_ */
