/*
 * FunctionMinimizer.hpp
 *
 *  Created on: Dec 11, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_FUNCTIONS_MINIMIZERS_FUNCTIONMINIMIZER_HPP_
#define MINIMIZATIONS_FUNCTIONS_MINIMIZERS_FUNCTIONMINIMIZER_HPP_

#include "BassoConfig.h"

struct FunctionMinimizer
{
	FunctionMinimizer(const unsigned int _N) :
		maxiterations(100),
		optimum(0),
		tempoptimum(_N),
		tempgradient(_N)
	{}

	//!> typedef for instance wrapped in shared pointer
	typedef boost::shared_ptr<FunctionMinimizer> ptr_t;

	typedef std::vector<double> array_type;

	//!> enumerates possible return codes of \sa CheckGradient().
	enum GradientStatus {
		gradient_success=0,	//!< gradient is (almost) zero
		gradient_continue=1,	//!< gradient is not yet zero
		error_noprogress=2, //!< gradient has not changed
		MAX_GradientStatus
	};

	//!> typedef for the function to check when to stop the iteration
	typedef boost::function<
			enum FunctionMinimizer::GradientStatus (const double) > check_function_t;

	/** Default check whether iteration should terminate as gradient's
	 * magnitude is below tolerance threshold \a _tol.
	 *
	 * @param _tol tolerance threshold
	 * @return status code
	 */
	virtual enum FunctionMinimizer::GradientStatus checkGradient(const double _tol) const = 0;

	typedef boost::function<double (
					const array_type &x)> function_evaluator_t;

	/** Setter for the external function to minimize.
	 *
	 * @param _fct function pointer
	 */
	void setFunctionEvaluator(const function_evaluator_t _fct) {
		function_evaluator = _fct;
	}

	typedef boost::function<array_type (
					const array_type &x)> gradient_evaluator_t;

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
	virtual void setMaxIterations(const unsigned int _maxiterations) = 0;

	/** Minimizes the specific functions with its gradients.
	 *
	 * @param _Tol desired error threshold
	 * @param _startvalue initial value and optimum on return
	 * @param _checkfunction check function in iteration
	 * @return number of iterations till threshold
	 */
	virtual const unsigned int minimize(
			const double _Tol,
			array_type &_startvalue,
			const check_function_t &_checkfunction
			) = 0;

	/** Getter for the value of the current optimum argument.
	 *
	 * @return optimum value
	 */
	const double getCurrentOptimumValue() const
	{ return optimum; }

	/** Getter for the current optimum argument.
	 *
	 * @return current optimum argument
	 */
	const array_type &getCurrentOptimum() const
	{ return tempoptimum; }

	/** Getter for the gradient at the current optimum argument.
	 *
	 * @return gradient at current optimum argument
	 */
	const array_type &getCurrentGradient() const
	{ return tempgradient; }

protected:
	unsigned int maxiterations;

	mutable double optimum;
	mutable std::vector<double> tempoptimum;
	mutable std::vector<double> tempgradient;

	//!> internal bound function to evaluate function itself
	function_evaluator_t function_evaluator;
	//!> internal bound function to evaluate gradient of function
	gradient_evaluator_t gradient_evaluator;
};



#endif /* MINIMIZATIONS_FUNCTIONS_MINIMIZERS_FUNCTIONMINIMIZER_HPP_ */
