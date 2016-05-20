/*
 * FunctionalMinimizer_inexactLinesearch_impl.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef FUNCTIONALMINIMIZER_INEXACTLINESEARCH_IMPL_HPP_
#define FUNCTIONALMINIMIZER_INEXACTLINESEARCH_IMPL_HPP_

#include "BassoConfig.h"

#include "Math/Helpers.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer_inexactLinesearch.hpp"

#include <boost/bind.hpp>

#include "Log/Logging.hpp"

template <class S>
const typename FunctionalMinimizer_inexactLinesearch<S>::Wolfe_indexset_t
FunctionalMinimizer_inexactLinesearch<S>::emptyset;

template <class S>
FunctionalMinimizer_inexactLinesearch<S>::FunctionalMinimizer_inexactLinesearch(
		const MinimizationFunctional<S> &_functional,
		FunctionMinimizer::ptr_t &_minimizer,
		S &_value) :
		FunctionalMinimizer<S>::FunctionalMinimizer(
				_functional, _minimizer, _value),
		Wolfe_indexset(FunctionalMinimizer_inexactLinesearch<S>::emptyset),
		constant_positivity(1e-4),
		constant_interpolation(1.)
{}

template <class S>
enum FunctionMinimizer::GradientStatus
FunctionalMinimizer_inexactLinesearch<S>::checkWolfeConditions(
		const double _startvalue,
		const array_type & _startgradient,
		const Wolfe_indexset_t &_Wolfe_indexset,
		const double _Tol) const
{
	const array_type currentiterate =
			minimizer->getCurrentOptimum();
	const double currentvalue =
			minimizer->getCurrentOptimumValue();
	const array_type currentgradient =
			minimizer->getCurrentGradient();

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
//		BOOST_LOG_TRIVIAL(trace)
//			<< "1. Positivity: " << component << " > "
//			<< constant_positivity;
//		conditions_fulfilled &=
//				component > constant_positivity;
//		// 2. still descent in the current gradient component
//		BOOST_LOG_TRIVIAL(trace)
//			<< "2. Still descent: " << currentgradient[*iter]
//			<< " < " << std::numeric_limits<double>::epsilon();
//		conditions_fulfilled &= currentgradient[*iter]
//						< std::numeric_limits<double>::epsilon();
//		// 3. stronger descent than linear interpolation
//		const double linearinterpolate = _startvalue
//				+ constant_interpolation * component *
//				_startgradient[*iter];
//		BOOST_LOG_TRIVIAL(trace)
//			<< "3. Stronger than linear: " << currentvalue
//			<< " < " << linearinterpolate;
//		conditions_fulfilled &= (currentvalue < linearinterpolate);
	if (conditions_fulfilled)
		return FunctionMinimizer::gradient_success;
	else
		return FunctionMinimizer::gradient_continue;
}

template <class S>
const unsigned int
FunctionalMinimizer_inexactLinesearch<S>::operator()(
		const unsigned int _N,
		const double _Tol,
		S &_startvalue) const
{
	std::vector<double> zeroposition(_N, 0.);
	const double functionzero =
			FunctionalMinimizer<S>::FunctionCaller(zeroposition);
	const std::vector<double> zerogradient =
			FunctionalMinimizer<S>::GradientCaller(zeroposition);

	// use Wolfe conditions to stop line search
	typedef boost::function<
			enum FunctionMinimizer::GradientStatus (const double) > check_function_t;
	check_function_t checkfunction =
			boost::bind(&FunctionalMinimizer_inexactLinesearch<S>::checkWolfeConditions,
					boost::cref(*this),
					functionzero,
					boost::cref(zerogradient),
					boost::cref(Wolfe_indexset),
					_1);

	array_type functionargument(_N, 0.);
	functional.convertInternalTypeToArrayType(_startvalue, functionargument);
	return minimizer->minimize(_Tol, functionargument, checkfunction);
}

#endif /* FUNCTIONALMINIMIZER_INEXACTLINESEARCH_IMPL_HPP_ */
