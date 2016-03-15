/*
 * FunctionalMinimizer_exactLinesearch.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef FUNCTIONALMINIMIZER_EXACT_LINESEARCH_HPP_
#define FUNCTIONALMINIMIZER_EXACT_LINESEARCH_HPP_

#include "BassoConfig.h"

#include "Minimizations/Functions/Minimizers/FunctionMinimizer.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizer.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Functions/Minimizers/Minimizer.hpp"

template <class S>
class FunctionalMinimizer_exactLinesearch : public FunctionalMinimizer<S>
{
	// bring some names of FunctionalMinimizer into our scope
	using FunctionalMinimizer<S>::functional;
	using FunctionalMinimizer<S>::minimizer;
	using FunctionalMinimizer<S>::value;

	using typename FunctionalMinimizer<S>::array_type;

public:
	//!> typedef for this instance wrapped in shared ptr
	typedef boost::shared_ptr< FunctionalMinimizer_exactLinesearch<S> > ptr_t;

	/** Constructor of functor FunctionalMinimizer_exactLinesearch.
	 *
	 * \note we use default values for Newton methods as suggested by
	 * [Nocedal, Wright, '99].
	 *
	 * @param _functional the function to minimize
	 * @param _minimizer the minimizer doing the minimization
	 * @param _value the external value as workspace for minimization
	 */
	FunctionalMinimizer_exactLinesearch(
			const MinimizationFunctional<S> &_functional,
			FunctionMinimizer::ptr_t &_minimizer,
			S &_value);

	//!> typedef for the function to check when to stop the iteration
	typedef boost::function<
			enum FunctionMinimizer::GradientStatus (const double) > check_function_t;

	/** Minimizes the given \a functional.
	 *
	 * @param _N dimension of the vector
	 * @param _Tol tolerance for minimization
	 * @param _startvalue initial value and minimizer in output
	 * @return required iterations
	 */
	const unsigned int operator()(
			const unsigned int _N,
			const double _Tol,
			S &_startvalue) const;
};

#include "FunctionalMinimizer_exactLinesearch_impl.hpp"


#endif /* FUNCTIONALMINIMIZER_EXACT_LINESEARCH_HPP_ */
