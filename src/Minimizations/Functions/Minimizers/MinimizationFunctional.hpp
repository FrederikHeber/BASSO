/*
 * MinimizationFunctional.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONFUNCTIONAL_HPP_
#define MINIMIZATIONFUNCTIONAL_HPP_

#include <vector>

/** This class defines the interface for a function
 * \f$ f: T \rightarrow \reel \f$ to be minimized.
 *
 * Here, T is some multi-dimensional object that must be convertible
 * into a array_type.
 *
 * This needs to be implemented for use with class \ref FunctionalMinimizer
 */
template <class T>
struct MinimizationFunctional
{
	//!> typedef for the internal array type for which conversion functions must be given
	typedef std::vector<double> array_type;

	/** Evaluates the function at \a _value.
	 *
	 * @param _value value where to evaluate function
	 * @return function value
	 */
	virtual double function(const T &_value) const = 0;

	/** Evaluates the function's gradient at \a _value.
	 *
	 * @param _value value where to evaluate function's gradient
	 * @return gradient direction
	 */
	virtual const T gradient(const T &_value) const = 0;

	/** We need to know how to convert the array type into a template type.
	 *
	 * Hence, this functions needs to be implemented.
	 *
	 * @param _x ref to array type
	 * @param _t value of the template type
	 */
	virtual void convertArrayTypeToInternalType(
			const array_type & _x,
			T &_t
			) const = 0;

	/** We need to know how to convert the template type to a array type.
	 *
	 * Hence, this functions needs to be implemented.
	 *
	 * @param _x ref to array type
	 * @param _t value of the template type
	 */
	virtual void convertInternalTypeToArrayType(
			const T &_t,
			array_type & _x
			) const = 0;
};

#endif /* MINIMIZATIONFUNCTIONAL_HPP_ */
