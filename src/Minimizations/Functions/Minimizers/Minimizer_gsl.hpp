/*
 * Minimizer_gsl.hpp
 *
 *  Created on: Nov 4, 2014
 *      Author: heber
 */

#ifndef MINIMIZER_GSL_HPP_
#define MINIMIZER_GSL_HPP_

#include "BassoConfig.h"

#include "Minimizer.hpp"

#include <boost/function.hpp>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>

#include <vector>

#include "Minimizations/Functions/Minimizers/FunctionMinimizer.hpp"

template <>
class Minimizer<gsl_vector> : public FunctionMinimizer
{
	typedef std::vector<double> array_type;

public:
	/** Constructor for class Minimizer.
	 *
	 * Allocates memory for the internally used minimization library.
	 *
	 * @param _N dimension, i.e. number of arguments of the function to minimize
	 */
	Minimizer(const unsigned int _N);

	/** Destructor for class Minimizer.
	 *
	 * Frees memory used by the internally used minimization library.
	 *
	 */
	virtual ~Minimizer();

	/** Default check whether iteration should terminate as gradient's
	 * magnitude is below tolerance threshold \a _tol.
	 *
	 * @param _tol tolerance threshold
	 * @return status code
	 */
	enum FunctionMinimizer::GradientStatus checkGradient(const double _tol) const;

	/** Setter for the external function to minimize.
	 *
	 * @param _fct function pointer
	 */
	void setFunctionEvaluator(const function_evaluator_t _fct) {
		function_evaluator = _fct;
	}

	/** Setter for the gradient of the external funciton to minimize.
	 *
	 * @param _fct gradient function pointer
	 */
	void setGradientEvaluator(const gradient_evaluator_t _fct) {
		gradient_evaluator = _fct;
	}

	/** Specifies a maximum number of iterations after which to stop.
	 *
	 * @param _maxiterations maximum number of iterations.
	 */
	void setMaxIterations(const unsigned int _maxiterations)
	{ maxiterations = _maxiterations; }

	//!> typedef for the function to check when to stop the iteration
	typedef boost::function<
			enum FunctionMinimizer::GradientStatus (const double) > check_function_t;

	/** Minimizes the specific functions with its gradients.
	 *
	 * @param _Tol desired error threshold
	 * @param _startvalue initial value and optimum on return
	 * @param _checkfunction check function in iteration
	 * @return number of iterations till threshold
	 */
	const unsigned int minimize(
			const double _Tol,
			array_type &_startvalue,
			const check_function_t &_checkfunction
			);

private:
	/** Converter for array type to the minimizer's internal type.
	 *
	 * @param _t array type instance
	 * @param _x internal type instance
	 */
	static void convertArrayTypeToInternalType(
			const array_type & _t,
			gsl_vector * const _x);

	/** Converter for the minimizer's internal type  into the array type.
	 *
	 * @param _x internal type instance
	 * @param _t array type instance
	 */
	static void convertInternalTypeToArrayType(
			const gsl_vector * const _x,
			array_type &_t);

	/** Wrapper for function_evaluator_t to translate array type into
	 * internal argument type for the specific minimizer.
	 *
	 * This needs to have the correct signature for the minimization
	 * library to be admissable as function pointer.
	 *
	 * @param x internal type representing the function's argument
	 * @param adata Minimizer instance to access function_evaluator_t
	 * @return function value
	 */
	static double FunctionCaller(
			const gsl_vector *x,
			void *adata);

	/** Wrapper for gradient_evaluator_t to translate array type into
	 * internal argument type for the specific minimizer.
	 *
	 * This needs to have the correct signature for the minimization
	 * library to be admissable as function pointer.
	 *
	 * @param x internal type representing the function's argument
	 * @param adata Minimizer instance to access gradient_evaluator_t
	 * @param g internal type representing the gradient on return
	 */
	static void GradientCaller(
			const gsl_vector *x,
			void *adata,
			gsl_vector *g);

	/** Wrapper for functiongradient_evaluator_t to translate array type
	 * into internal argument type for the specific minimizer.
	 *
	 * This needs to have the correct signature for the minimization
	 * library to be admissable as function pointer.
	 *
	 * @param x internal type representing the function's argument
	 * @param adata Minimizer instance to access functiongradient_evaluator_t
	 * @param f function value on return
	 * @param g internal type representing the gradient on return
	 */
	static void FunctionGradientCaller(
			const gsl_vector *x,
			void *adata,
			double *f,
			gsl_vector *g);

private:
	const gsl_multimin_fdfminimizer_type *minimizer;
	gsl_multimin_fdfminimizer *s;
	gsl_multimin_function_fdf my_func;
	const unsigned int N;
};

#endif /* MINIMIZER_GSL_HPP_ */
