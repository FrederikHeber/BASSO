/*
 * FunctionalMinimizer.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef FUNCTIONALMINIMIZER_HPP_
#define FUNCTIONALMINIMIZER_HPP_

#include "BassoConfig.h"

#include "Minimizations/Functions/Minimizers/Minimization_common.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Functions/Minimizers/Minimizer.hpp"

#include <vector>

#include <gsl/gsl_multimin.h>

//!> index set for Wolfe conditions
typedef std::vector<unsigned int> Wolfe_indexset_t;

/** This class is a wrapper to library minimization routines.
 *
 * To use it you need to derive MinimizationFunctional<S> for your
 * specific type T and implement the virtual functions.
 *
 * Second, choose a minimizer specialization Minimizer<T> such as
 * Minimizer<gsl_vector> using GSL routines.
 *
 * \code
 * #include "FunctionalMinimizer.hpp"
 * #include "MinimizationFunctional.hpp"
 * #include "Minimizer.hpp"
 *
 * class MyFunctional :
 *   public MinimizationFunctional<MyType>
 * {
 * 		virtual double operator()(const MyType &_value) const
 * 		{
 * 		  // ...
 * 		}
 * 		// ... and also the other functions ...
 * };
 *
 * MyFunctional functional; // the instantiated derived functional
 * unsigned int dim = 4;
 * Minimizer<gsl_vector> minimizer(dim);	// library minimizer itself
 * // a temporary value for the minizer to work on with parameter dim to init
 * MyType temp(dim);
 * //initializing the minimizer with the functional and temporary
 * FunctionalMinimizer<MyType> fmin(functional, minimizer, temp);
 * \endcode
 *
 * Then, you may finally call the minimizer to do its work:
 * \code
 * MyType startvalue = MyTypeZeroValue(dim); // zero as start value in MyType
 * MyType min = fmin(dim, 1e-6, startvalue);
 * \endcode
 *
 * The resulting minimizer of the function in MyFunctional is obtained
 * by minimizer::getCurrentOptimum().
 *
 * \section wolfe-conditions Wolfe conditions
 *
 * The function minimization can also be performed inexact. Then Wolfe
 * conditions for a specific subset of indices (_value is a vector, right?)
 * and if they are fulfilled iteration is discontinued without any regard
 * to magnitude of gradient or else.
 *
 * For this you have to use setinexactLinesearch() and give it a suitable
 * Wolfe_indexset_t with the indices that are truly descent directions and for
 * which the three Wolfe conditions - positivity, still descent and stronger
 * descent than the linear interpolation -- must hold.
 */
template <class S, class T>
class FunctionalMinimizer
{
	typedef typename MinimizationFunctional<S>::array_type array_type;
public:
	/** Constructor of functor FunctionalMinimizer.
	 *
	 * \note we use default values for Newton methods as suggested by
	 * [Nocedal, Wright, '99].
	 *
	 * @param _functional the function to minimize
	 * @param _value the external value as workspace for minimization
	 * @param _constant_positivity Wolfe constant for positivity
	 * @param _constant_interpolation Wolfe constant for stronger than linear
	 */
	FunctionalMinimizer(
			const MinimizationFunctional<S> &_functional,
			Minimizer<T> &_minimizer,
			S &_value,
			const double _constant_positivity = 1e-4,
			const double _constant_interpolation = 1.);

	//!> typedef for the function to check when to stop the iteration
	typedef boost::function<
			enum Minimization::GradientStatus (const double) > check_function_t;

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

	/** Minimizes the given \a functional with an inexact
	 * linesearch fulfilling Wolfe conditions on a (sub)set of indices.
	 *
	 * @param _N dimension of the vector
	 * @param _Tol tolerance for minimization
	 * @param _Wolfe_indexset (sub)set of indices (i.e. those directions that
	 *        are truly descent directions)
	 * @param _startvalue initial value and minimizer in output
	 * @return required iterations
	 */
	const unsigned int operator()(
			const unsigned int _N,
			const double _Tol,
			const Wolfe_indexset_t &_Wolfe_indexset,
			S &_startvalue) const;

	/** Setter for the maximum number of iterations.
	 *
	 * @param _maxiterations maximum number of iterations
	 */
	void setMaxIterations(const unsigned int _maxiterations)
	{ minimizer.setMaxIterations(_maxiterations); }

	/** Setter for inexact line searches.
	 *
	 * We stop the iteration as soon as the Wolfe conditions are fulfilled
	 * for the given (sub)set of indices \a _Wolfe_indexset (value is a
	 * vector, right?).
	 *
	 * @param _inexactLinesearch true - do inexact line search with Wolfe
	 *        conditions
	 * @param _Wolfe_indexset subset of indices (i.e. the directions that
	 *        are descent directions)
	 */
	void setinexactLinesearch(
			const bool _inexactLinesearch,
			const Wolfe_indexset_t &_Wolfe_indexset)
	{
		inexactLinesearch = _inexactLinesearch;
		Wolfe_indexset = _Wolfe_indexset;
	}

	//!> instance of the minimization functional
	const MinimizationFunctional<S> &functional;
	//!> instance of the minimizer
	Minimizer<T> &minimizer;
	S &value;

private:

	/** Evaluates the Wolfe conditions for a certain step width \a _x with
	 * and (sub)set of indices \a _Wolfe_indexset.
	 *
	 * @param _startvalue function value with step width zero
	 * @param _startgradient gradient with step width zero
	 * @param _Wolfe_indexset (sub)set of indices (i.e. the descent directions)
	 * @param _Tol tolerance *unused*
	 * @return GradientStatus
	 */
	enum Minimization::GradientStatus checkWolfeConditions(
			const double _startvalue,
			const array_type & _startgradient,
			const Wolfe_indexset_t &_Wolfe_indexset,
			const double _Tol) const;

	typedef boost::function<double (
					const array_type &x)> function_evaluator_t;

	/** Wrapper for the call to MinimizationFunctional::function().
	 *
	 * @param x argument
	 * @return evaluated function
	 */
	double FunctionCaller(
			const array_type &x) const;

	typedef boost::function<array_type (
					const array_type &x)> gradient_evaluator_t;

	/** Wrapper for the call to MinimizationFunctional::gradient().
	 *
	 * @param x argument
	 * @return evaluated gradient
	 */
	array_type GradientCaller(
			const array_type &x) const;

private:
	// inexact line search values

	//!> whether to perform inexact line search with Wolfe conditions
	bool inexactLinesearch;
	//!> index set for which to check the Wolfe conditions
	Wolfe_indexset_t Wolfe_indexset;
	//!> constant above which the step width must always lie
	const double constant_positivity;
	//!> constant to scale the linear interpolation
	const double constant_interpolation;
};

#include "FunctionalMinimizer_impl.hpp"


#endif /* FUNCTIONALMINIMIZER_HPP_ */
