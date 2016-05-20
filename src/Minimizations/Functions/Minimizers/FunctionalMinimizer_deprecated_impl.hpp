/*
 * FunctionalMinimizer_deprecated_impl.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef FUNCTIONALMINIMIZER_DEPRECATED_IMPL_HPP_
#define FUNCTIONALMINIMIZER_DEPRECATED_IMPL_HPP_

#include "BassoConfig.h"

#include "Math/Helpers.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer_deprecated.hpp"

#include <boost/bind.hpp>

#include <boost/log/trivial.hpp>

template <class S, class T>
FunctionalMinimizer_deprecated<S,T>::FunctionalMinimizer_deprecated(
		const MinimizationFunctional<S> &_functional,
		Minimizer<T> &_minimizer,
		S &_value,
		const double _constant_positivity,
		const double _constant_interpolation) :
	functional(_functional),
	minimizer(_minimizer),
	value(_value),
	inexactLinesearch(false),
	constant_positivity(_constant_positivity),
	constant_interpolation(_constant_interpolation)
{
	minimizer.setFunctionEvaluator(
			boost::bind(&FunctionalMinimizer_deprecated<S,T>::FunctionCaller,
					boost::ref(*this),
					_1));
	minimizer.setGradientEvaluator(
			boost::bind(&FunctionalMinimizer_deprecated<S,T>::GradientCaller,
					boost::ref(*this),
					_1));
}

template <class S, class T>
enum FunctionMinimizer::GradientStatus
FunctionalMinimizer_deprecated<S,T>::checkWolfeConditions(
		const double _startvalue,
		const array_type & _startgradient,
		const Wolfe_indexset_t &_Wolfe_indexset,
		const double _Tol) const
{
	const array_type currentiterate =
			minimizer.getCurrentOptimum();
	const double currentvalue =
			minimizer.getCurrentOptimumValue();
	const array_type currentgradient =
			minimizer.getCurrentGradient();

	bool conditions_fulfilled = true;
	assert( !_Wolfe_indexset.empty() );

	double linearinterpolate = _startvalue;
	double interpolatedgradient = 0.;
	double realgradient = 0.;
	Wolfe_indexset_t indices(currentiterate.size());
	std::generate( indices.begin(), indices.end(), Helpers::unique_number());
	for (Wolfe_indexset_t::const_iterator iter = _Wolfe_indexset.begin();
			(iter != _Wolfe_indexset.end()) && conditions_fulfilled;
			++iter) {
		const double componentiterate = currentiterate[*iter];
		const double componentgradient = currentgradient[*iter];

		linearinterpolate += -constant_positivity * componentiterate *
				::pow(_startgradient[*iter],2);

		interpolatedgradient += constant_interpolation
				* ::pow(_startgradient[*iter],2);
		realgradient += componentgradient * _startgradient[*iter];
	}
	// 1. sufficient decrease
	LOG(debug, "1. sufficient decrease: " << currentvalue << " <= " << linearinterpolate);
	conditions_fulfilled &=
			(currentvalue <= linearinterpolate);

	// 2. curvature condition
	LOG(debug, "2. curvature condition: " << realgradient << " <= " << interpolatedgradient);
	conditions_fulfilled &=
			realgradient <= interpolatedgradient;

//		// 1. positivity of step width component
//		LOG(trace, "1. Positivity: " << component << " > " << constant_positivity);
//		conditions_fulfilled &=
//				component > constant_positivity;
//		// 2. still descent in the current gradient component
//		LOG(trace, "2. Still descent: " << currentgradient[*iter]
//			<< " < " << std::numeric_limits<double>::epsilon());
//		conditions_fulfilled &= currentgradient[*iter]
//						< std::numeric_limits<double>::epsilon();
//		// 3. stronger descent than linear interpolation
//		const double linearinterpolate = _startvalue
//				+ constant_interpolation * component *
//				_startgradient[*iter];
//		LOG(trace, "3. Stronger than linear: " << currentvalue
//			<< " < " << linearinterpolate);
//		conditions_fulfilled &= (currentvalue < linearinterpolate);
	if (conditions_fulfilled)
		return FunctionMinimizer::gradient_success;
	else
		return FunctionMinimizer::gradient_continue;
}

template <class S, class T>
const unsigned int
FunctionalMinimizer_deprecated<S,T>::operator()(
		const unsigned int _N,
		const double _Tol,
		S &_startvalue) const
{
	// use the default check to stop line search
	check_function_t checkfunction =
			boost::bind(&Minimizer<T>::checkGradient,
					boost::cref(minimizer), _1);

	array_type functionargument(_N, 0.);
	functional.convertInternalTypeToArrayType(_startvalue, functionargument);
	const unsigned int iterations =
			minimizer.minimize(_Tol, functionargument, checkfunction);
	functional.convertArrayTypeToInternalType(functionargument, _startvalue);
	return iterations;
}

template <class S, class T>
const unsigned int
FunctionalMinimizer_deprecated<S,T>::operator()(
		const unsigned int _N,
		const double _Tol,
		const Wolfe_indexset_t &_Wolfe_indexset,
		S &_startvalue) const
{
	std::vector<double> zeroposition(_N, 0.);
	const double functionzero = FunctionCaller(
			zeroposition);
	const std::vector<double> zerogradient =
			GradientCaller(zeroposition);

	// use Wolfe conditions to stop line search
	check_function_t checkfunction =
			boost::bind(&FunctionalMinimizer_deprecated<S,T>::checkWolfeConditions,
					boost::cref(*this),
					functionzero,
					boost::cref(zerogradient),
					boost::cref(_Wolfe_indexset),
					_1);

	array_type functionargument(_N, 0.);
	functional.convertInternalTypeToArrayType(_startvalue, functionargument);
	return minimizer.minimize(_Tol, functionargument, checkfunction);
}

template <class S, class T>
double FunctionalMinimizer_deprecated<S,T>::FunctionCaller(
		const array_type &x) const
{
	functional.convertArrayTypeToInternalType(x, value);
	return functional.function(value);
}

template <class S, class T>
typename FunctionalMinimizer_deprecated<S,T>::array_type
FunctionalMinimizer_deprecated<S,T>::GradientCaller(
		const array_type &x) const
{
	functional.convertArrayTypeToInternalType(x, value);
	const S tempgradient = functional.gradient(value);
	array_type returngradient(x.size(), 0.);
	functional.convertInternalTypeToArrayType(tempgradient, returngradient);
	return returngradient;
}


#endif /* FUNCTIONALMINIMIZER_DEPRECATED_IMPL_HPP_ */
