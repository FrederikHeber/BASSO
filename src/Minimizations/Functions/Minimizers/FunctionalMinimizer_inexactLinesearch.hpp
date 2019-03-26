/*
 * FunctionalMinimizer_inexactLinesearch.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef FUNCTIONALMINIMIZER_INEXACT_LINESEARCH_HPP_
#define FUNCTIONALMINIMIZER_INEXACT_LINESEARCH_HPP_

#include "BassoConfig.h"

#include "Minimizations/Functions/Minimizers/FunctionalMinimizer.hpp"
#include "Minimizations/Functions/Minimizers/FunctionMinimizer.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Functions/Minimizers/Minimizer.hpp"

#include <vector>

class FunctionalMinimizerFactory;

template <class S>
class FunctionalMinimizer_inexactLinesearch : public FunctionalMinimizer<S>
{
	// bring some names of FunctionalMinimizer into our scope
	using FunctionalMinimizer<S>::functional;
	using FunctionalMinimizer<S>::minimizer;
	using FunctionalMinimizer<S>::value;

	using typename FunctionalMinimizer<S>::array_type;

public:
	//!> typedef for this instance wrapped in shared ptr
	typedef boost::shared_ptr< FunctionalMinimizer_inexactLinesearch<S> > ptr_t;

	/** Constructor of functor FunctionalMinimizer_inexactLinesearch.
	 *
	 * \note we use default values for Newton methods as suggested by
	 * [Nocedal, Wright, '99].
	 *
	 * @param _functional the function to minimize
	 * @param _minimizer the minimizer doing the minimization
	 * @param _value the external value as workspace for minimization
	 */
	FunctionalMinimizer_inexactLinesearch(
			const MinimizationFunctional<S> &_functional,
			FunctionMinimizer::ptr_t &_minimizer,
			S &_value);

	virtual ~FunctionalMinimizer_inexactLinesearch() {}

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

private:
	//! grant factory access to private setter
	friend class FunctionalMinimizerFactory;

	//!> index set for Wolfe conditions
	typedef std::vector<unsigned int> Wolfe_indexset_t;

	//!> define a static empty Wolfe index set as default
	static const Wolfe_indexset_t emptyset;

	/** Setter for inexact line searches.
	 *
	 * We stop the iteration as soon as the Wolfe conditions are fulfilled
	 * for the given (sub)set of indices \a _Wolfe_indexset (value is a
	 * vector, right?).
	 *
	 * @param _Wolfe_indexset subset of indices (i.e. the directions that
	 *        are descent directions)
	 * @param _constant_positivity Wolfe constant for positivity
	 * @param _constant_interpolation Wolfe constant for stronger than linear
	 */
	void setInexactLinesearchParameters(
			const Wolfe_indexset_t &_Wolfe_indexset,
			const double _constant_positivity = 1e-4,
			const double _constant_interpolation = 1.)
	{
		const_cast<Wolfe_indexset_t &>(Wolfe_indexset) = _Wolfe_indexset;
		const_cast<double &>(constant_positivity) = _constant_positivity;
		const_cast<double &>(constant_interpolation) = _constant_interpolation;
	}

	/** Evaluates the Wolfe conditions for a certain step width \a _x with
	 * and (sub)set of indices \a _Wolfe_indexset.
	 *
	 * @param _startvalue function value with step width zero
	 * @param _startgradient gradient with step width zero
	 * @param _Wolfe_indexset (sub)set of indices (i.e. the descent directions)
	 * @param _Tol tolerance *unused*
	 * @return GradientStatus
	 */
	enum FunctionMinimizer::GradientStatus checkWolfeConditions(
			const double _startvalue,
			const array_type & _startgradient,
			const Wolfe_indexset_t &_Wolfe_indexset,
			const double _Tol) const;

private:

	//!> index set for which to check the Wolfe conditions
	const Wolfe_indexset_t &Wolfe_indexset;
	//!> constant above which the step width must always lie
	const double constant_positivity;
	//!> constant to scale the linear interpolation
	const double constant_interpolation;
};

#include "FunctionalMinimizer_inexactLinesearch_impl.hpp"


#endif /* FUNCTIONALMINIMIZER_INEXACT_LINESEARCH_HPP_ */
