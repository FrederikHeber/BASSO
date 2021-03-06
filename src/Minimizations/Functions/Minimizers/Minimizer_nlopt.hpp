/*
 * Minimizer_nlopt.hpp
 *
 *  Created on: Nov 6, 2014
 *      Author: heber
 */

#ifndef MINIMIZER_NLOPT_HPP_
#define MINIMIZER_NLOPT_HPP_

#include "BassoConfig.h"

#include "Minimizer.hpp"

#include <boost/function.hpp>

#include <vector>

#include "Minimizations/Functions/Minimizers/FunctionMinimizer.hpp"

#ifdef NLOPT_FOUND
#include <nlopt.hpp>
#endif /* NLOPT_FOUND */

typedef std::vector<double> NLopt_vector;

template <>
class Minimizer<NLopt_vector> : public FunctionMinimizer
{
public:
	/** Constructor for class Minimizer.
	 *
	 * Allocates memory for the internally used minimization library.
	 *
	 * @param _N dimension, i.e. number of arguments of the function to minimize
	 */
	Minimizer(const unsigned int _N);

	virtual ~Minimizer() {}

	typedef boost::function<double (
			const std::vector<double> &,	/* x */
			std::vector<double> &, 	/* grad */
			void *my_func_data
			)> functiongradient_evaluator_t;

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

	void setPositivityBoundIndices(
			const std::vector<unsigned int> &_indices)
	{ PositivityBoundIndices = _indices; }

	void setConstantPositivity(
			const double _constant)
	{ constant_positivity = _constant; }

	/** Default check whether iteration should terminate as gradient's
	 * magnitude is below tolerance threshold \a _tol.
	 *
	 * @param _tol tolerance threshold
	 * @return status code
	 */
	enum FunctionMinimizer::GradientStatus checkGradient(const double _tol) const;

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
			NLopt_vector * const _x);

	/** Converter for the minimizer's internal type  into the array type.
	 *
	 * @param _x internal type instance
	 * @param _t array type instance
	 */
	static void convertInternalTypeToArrayType(
			const NLopt_vector * const _x,
			array_type &_t);

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
	static double FunctionGradientCaller(
			const std::vector<double> &,	/* x */
			std::vector<double> &, 	/* grad */
			void *my_func_data
			);

private:
#ifdef NLOPT_FOUND
	nlopt::opt opt;
#endif /* NLOPT_FOUND */
	const unsigned int N;
	unsigned int iter;

	std::vector<unsigned int> PositivityBoundIndices;
	double constant_positivity;

	mutable check_function_t checkfunction;
	mutable double tolerance;
};

#endif /* MINIMIZER_NLOPT_HPP_ */
