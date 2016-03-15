/*
 * FunctionalMinimizer_impl.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef FUNCTIONALMINIMIZER_IMPL_HPP_
#define FUNCTIONALMINIMIZER_IMPL_HPP_

#include "BassoConfig.h"

#include "Math/Helpers.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer.hpp"

#include <boost/bind.hpp>

#include <boost/log/trivial.hpp>

template <class S>
FunctionalMinimizer<S>::FunctionalMinimizer(
		const MinimizationFunctional<S> &_functional,
		FunctionMinimizer::ptr_t &_minimizer,
		S &_value) :
	functional(_functional),
	minimizer(_minimizer),
	value(_value)
{
	minimizer->setFunctionEvaluator(
			boost::bind(&FunctionalMinimizer<S>::FunctionCaller,
					boost::ref(*this),
					_1));
	minimizer->setGradientEvaluator(
			boost::bind(&FunctionalMinimizer<S>::GradientCaller,
					boost::ref(*this),
					_1));
}

template <class S>
double FunctionalMinimizer<S>::FunctionCaller(
		const array_type &x) const
{
	functional.convertArrayTypeToInternalType(x, value);
	return functional.function(value);
}

template <class S>
typename FunctionalMinimizer<S>::array_type
FunctionalMinimizer<S>::GradientCaller(
		const array_type &x) const
{
	functional.convertArrayTypeToInternalType(x, value);
	const S tempgradient = functional.gradient(value);
	array_type returngradient(x.size(), 0.);
	functional.convertInternalTypeToArrayType(tempgradient, returngradient);
	return returngradient;
}


#endif /* FUNCTIONALMINIMIZER_IMPL_HPP_ */
