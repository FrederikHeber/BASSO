/*
 * FunctionMinimizer.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef FUNCTIONMINIMIZER_HPP_
#define FUNCTIONMINIMIZER_HPP_

#include "BassoConfig.h"

#include "Minimizations/Functions/MinimizationFunctional.hpp"

/** This class is a wrapper to the gsl minimization routines.
 *
 * To use it you need to derive MinimizationFunctional<T> for your
 * specific type T and implement the virtual functions.
 *
 * \code
 * #include "FunctionMinimizer.hpp"
 * #include "MinimizationFunctional.hpp"
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
 * // a temporary value for the minizer to work on with parameter dim to init
 * MyType temp(dim);
 * //initializing the minimizer with the functional and temporary
 * FunctionMinimizer<MyType> minimizer(functional, temp);
 * \endcode
 *
 * Finally, do not forget to instantiate template functions below via the
 * defined macro in your code, e.g.
 * \code
 * CONSTRUCT_FUNCTIONMINIMIZER(Eigen::VectorXd)
 * \endcode
 *
 * Then, you may finally call the minimizer to do its work:
 * \code
 * MyType startvalue = MyTypeZeroValue(dim); // zero as start value in MyType
 * MyType min = minimizer(dim, 1e-6, startvalue);
 * \endcode
 *
 * The resulting minimizer of the function in MyFunctional is in \a min.
 */
template <class T>
class FunctionMinimizer
{
public:
	FunctionMinimizer(
			const MinimizationFunctional<T> &_functional,
			T &_value) :
		functional(_functional),
		value(_value)
	{}
	~FunctionMinimizer() {}

	/** Performs the minimization on the given \a functional.
	 *
	 * @param _N dimension of the vector
	 * @param _Tol tolerance for minimization
	 * @param _startvalue initial value
	 * @return minimizer
	 */
	const T operator()(
			const unsigned int _N,
			const double _Tol,
			const T &_startvalue);

	//!> instance of the minimization functional
	const MinimizationFunctional<T> &functional;
	//!> internal (initialized!) value to work on
	T &value;
};

template <class T>
double
FunctionMinimizer_FunctionCaller(
		const gsl_vector *x,
		void *adata)
{
	// obtain minimizer object
	struct FunctionMinimizer<T> *minimizer =
			static_cast<FunctionMinimizer<T> *>(adata);
	// convert given value to something interpretable
	minimizer->functional.convertToInternalType(
			minimizer->value,
			x);
	// evaluate the function
	const double tmp =
			minimizer->functional(minimizer->value);
	// return the value
	return tmp;
}

template <class T>
void
FunctionMinimizer_GradientCaller(
		const gsl_vector *x,
		void *adata,
		gsl_vector *g)
{
	// obtain minimizer object
	struct FunctionMinimizer<T> *minimizer =
			static_cast<FunctionMinimizer<T> *>(adata);
	// convert given value to something interpretable
	minimizer->functional.convertToInternalType(
			minimizer->value,
			x);
	// evaluate the gradient
	const T tmpgradient =
			minimizer->functional.gradient(minimizer->value);
	// convert gradient and store in g
	minimizer->functional.convertFromInternalType(tmpgradient, g);
}

template <class T>
void
FunctionMinimizer_FunctionGradientCaller(
		const gsl_vector *x,
		void *adata,
		double *f,
		gsl_vector *g)
{
	// obtain minimizer object
	struct FunctionMinimizer<T> *minimizer =
			static_cast<FunctionMinimizer<T> *>(adata);
	// convert given value to something interpretable
	minimizer->functional.convertToInternalType(minimizer->value, x);
	// evaluate the function and store in f
	*f = minimizer->functional(minimizer->value);
	// evaluate the gradient
	const T tmpgradient =
			minimizer->functional.gradient(minimizer->value);
	// convert gradient and store in g
	minimizer->functional.convertFromInternalType(tmpgradient, g);
}

//!> use this macro to make sure wrapper functions for gsl are instantiated
#define CONSTRUCT_FUNCTIONMINIMIZER(type) \
		template double FunctionMinimizer_FunctionCaller<type>(const gsl_vector *x,void *adata); \
		template void FunctionMinimizer_GradientCaller<type>(const gsl_vector *x,void *adata,gsl_vector *g); \
		template void FunctionMinimizer_FunctionGradientCaller<type>(const gsl_vector *x,void *adata,double *f,gsl_vector *g);

#include "Minimizations/Functions/FunctionMinimizer_impl.hpp"


#endif /* FUNCTIONMINIMIZER_HPP_ */
