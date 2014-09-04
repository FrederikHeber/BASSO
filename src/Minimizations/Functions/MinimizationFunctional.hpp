/*
 * MinimizationFunctional.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONFUNCTIONAL_HPP_
#define MINIMIZATIONFUNCTIONAL_HPP_

#include <gsl/gsl_vector.h>

/** This class defines the interface for a function
 * \f$ f: T \rightarrow \reel \f$ to be minimized.
 *
 * Here, T is some multi-dimensional object that must be convertible
 * into a gsl_vector.
 *
 * This needs to be implemented for use with class \ref FunctionMinimizer
 */
template <class T>
struct MinimizationFunctional
{
	/** Evaluates the function at \a _value.
	 *
	 * @param _value value where to evaluate function
	 * @return function value
	 */
	virtual double operator()(const T &_value) const = 0;

	/** Evaluates the function's gradient at \a _value.
	 *
	 * @param _value value where to evaluate function's gradient
	 * @return gradient direction
	 */
	virtual const T gradient(const T &_value) const = 0;

	/** We need to know how to convert a gsl_vector to the internal type.
	 *
	 * Hence, this functions needs to be implemented.
	 *
	 * @param _t value of the internal type
	 * @param x ptr to gsl_vector (which gsl minimization uses)
	 */
	virtual void convertToInternalType(
			T &_t,
			const gsl_vector * const x) const = 0;

	/** We need to know how to convert the internal type to a gsl_vector.
	 *
	 * Hence, this functions needs to be implemented.
	 *
	 * @param _t value of the internal type
	 * @param x ptr to gsl_vector (which gsl minimization uses)
	 */
	virtual void convertFromInternalType(
			const T &_t,
			gsl_vector * x) const = 0;
};

#endif /* MINIMIZATIONFUNCTIONAL_HPP_ */
