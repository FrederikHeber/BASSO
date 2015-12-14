/*
 * FunctionalMinimizer.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef FUNCTIONALMINIMIZER_HPP_
#define FUNCTIONALMINIMIZER_HPP_

#include "BassoConfig.h"

#include "Minimizations/Functions/Minimizers/FunctionMinimizer.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/Functions/Minimizers/Minimizer.hpp"

#include <vector>

//!> index set for Wolfe conditions
typedef std::vector<unsigned int> Wolfe_indexset_t;

/** This class is an interface to library minimization routines.
 *
 * To use it you need to derive MinimizationFunctional<S> for your
 * specific type T and implement the virtual functions.
 *
 * Second, choose a FunctionalMinimizer specialization by calling upon
 * FunctionalMinimizerFactory.
 *
 * \code
 * #include "FunctionalMinimizer.hpp"
 * #include "FunctionalMinimizerFactory"
 *
 *
 * MyFunctional functional; // the instantiated derived functional
 * unsigned int dim = 4;
 * std::vector<double> tmin(dim);
 * // instantiate the minimizer via the factory, handing over functional
 * // dim(ensions) and a place to store the minimum value tmin
 * FunctionalMinimizer< std::vector<double> >::ptr_t functionminimizer;
 *		functionminimizer =
 *				FunctionalMinimizerFactory::create< std::vector<double> >(
 *					dim,
 *					functional,
 *					tmin);
 * functionminimizer->setMaxIterations(MaxInnerIterations);
 * \endcode
 *
 * Then, you may finally call the minimizer to do its work:
 * \code
 * // call the minimizer
 * const double TolFun = 1-12; // tolerance value
 * unsigned int inner_iterations = (*functionminimizer)(
 *			dim, TolFun, tmin);
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
 * For this you have to use a different FunctionalMinimizerFactory::create()
 * function that also accepts parameters for specifying these Wolfe
 * conditions.
 */
template <class S>
class FunctionalMinimizer
{
protected:
	typedef typename MinimizationFunctional<S>::array_type array_type;
public:
	//!> typedef for this instance wrapped in shared ptr
	typedef boost::shared_ptr< FunctionalMinimizer<S> > ptr_t;

	/** Constructor of functor FunctionalMinimizer.
	 *
	 * \note we use default values for Newton methods as suggested by
	 * [Nocedal, Wright, '99].
	 *
	 * @param _functional the function to minimize
	 * @param _minimizer the minimizer doing the minimization
	 * @param _value the external value as workspace for minimization
	 */
	FunctionalMinimizer(
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
	virtual const unsigned int operator()(
			const unsigned int _N,
			const double _Tol,
			S &_startvalue) const = 0;

	/** Setter for the maximum number of iterations.
	 *
	 * @param _maxiterations maximum number of iterations
	 */
	void setMaxIterations(const unsigned int _maxiterations)
	{ minimizer->setMaxIterations(_maxiterations); }

	/** Getter for the value of the current optimum argument.
	 *
	 * @return optimum value
	 */
	const double getCurrentOptimumValue() const
	{ return minimizer->getCurrentOptimumValue(); }

	/** Getter for the current optimum argument.
	 *
	 * @return current optimum argument
	 */
	const FunctionMinimizer::array_type &getCurrentOptimum() const
	{ return minimizer->getCurrentOptimum(); }

	/** Getter for the gradient at the current optimum argument.
	 *
	 * @return gradient at current optimum argument
	 */
	const FunctionMinimizer::array_type &getCurrentGradient() const
	{ return minimizer->getCurrentGradient(); }

protected:

	//!> instance of the minimization functional
	const MinimizationFunctional<S> &functional;
	//!> instance of the minimizer
	FunctionMinimizer::ptr_t minimizer;
	S &value;

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
};

#include "FunctionalMinimizer_impl.hpp"


#endif /* FUNCTIONALMINIMIZER_HPP_ */
