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

#include <vector>

#include <gsl/gsl_multimin.h>

//!> index set for Wolfe conditions
typedef std::vector<unsigned int> Wolfe_indexset_t;

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
template <class T>
class FunctionMinimizer
{
	typedef typename MinimizationFunctional<T>::array_type array_type;
public:
	/** Constructor of functor FunctionMinimizer.
	 *
	 * \note we use default values for Newton methods as suggested by
	 * [Nocedal, Wright, '99].
	 *
	 * @param _functional the function to minimize
	 * @param _value the external value as workspace for minimization
	 * @param _constant_positivity Wolfe constant for positivity
	 * @param _constant_interpolation Wolfe constant for stronger than linear
	 */
	FunctionMinimizer(
			const MinimizationFunctional<T> &_functional,
			T &_value,
			const double _constant_positivity = 1e-4,
			const double _constant_interpolation = 1.) :
		functional(_functional),
		value(_value),
		maxiterations(100),
		inexactLinesearch(false),
		constant_positivity(_constant_positivity),
		constant_interpolation(_constant_interpolation)
	{}
	~FunctionMinimizer() {}

	/** Performs the minimization on the given \a functional.
	 *
	 * @param _N dimension of the vector
	 * @param _Tol tolerance for minimization
	 * @param _startvalue initial value
	 * @return required iterations
	 */
	const unsigned int operator()(
			const unsigned int _N,
			const double _Tol,
			T &_startvalue);

	/** Performs the minimization on the given \a functional with an inexact
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
			T &_startvalue);

	/** Setter for the maximum number of iterations.
	 *
	 * @param _maxiterations maximum number of iterations
	 */
	void setMaxIterations(const unsigned int _maxiterations)
	{ maxiterations = _maxiterations; }

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
	const MinimizationFunctional<T> &functional;
	//!> internal (initialized!) value to work on
	T &value;

private:
	/** Evaluates the Wolfe conditions for a certain step width \a _x with
	 * and (sub)set of indices \a _Wolfe_indexset.
	 *
	 * @param _startvalue function value with step width zero
	 * @param _startgradient gradient with step width zero
	 * @param _Wolfe_indexset (sub)set of indices (i.e. the descent directions)
	 * @param _s minimization structure with step width, function and gradient
	 *        values
	 * @return GSL_SUCCESS (stop) or GSL_CONTINUE (continue)
	 */
	int checkWolfeConditions(
			const double _startvalue,
			const gsl_vector * const _startgradient,
			const Wolfe_indexset_t &_Wolfe_indexset,
			const gsl_multimin_fdfminimizer * const _s) const;

	//!> typedef for the function to check when to stop the iteration
	typedef boost::function< int (
			const gsl_multimin_fdfminimizer * const,
			const double) > check_function_t;

	/** Helper function that contains the actula minimization procedure but
	 * is generalized to allow for given function to check when to stop the
	 * iteraton.
	 *
	 * @param _N dimensionality of function
	 * @param _Tol tolerance when to stop (for exact line search)
	 * @param _startvalue initial value for minimization
	 * @param _checkfunction bound function to check when to stop
	 * @return required iterations
	 */
	const unsigned int performMinimization(
			const unsigned int _N,
			const double _Tol,
			T &_startvalue,
			const check_function_t &_checkfunction
			);

public:
	/** We need to know how to convert the array_type to the internal
	 * type of the minimization library, here gsl_vector.
	 *
	 * Hence, this functions needs to be implemented.
	 *
	 * @param _t value of the internal type
	 * @param x ptr to gsl_vector (which gsl minimization uses)
	 */
	void convertToInternalType(
			const array_type & _t,
			gsl_vector * const _x) const;

	/** We need to know how to convert the internal type of the
	 * minimization to an array_type.
	 *
	 * Hence, this functions needs to be implemented.
	 *
	 * @param _t value of the internal type
	 * @param x ptr to gsl_vector (which gsl minimization uses)
	 */
	void convertFromInternalType(
			const gsl_vector * const _x,
			array_type &_t) const;

private:
	//!> set maximum number of iterations
	unsigned int maxiterations;

	// inexact line search values

	//!> whether to perform inexact line search with Wolfe conditions
	bool inexactLinesearch;
	//!> index set for which to check the Wolfe conditions
	Wolfe_indexset_t Wolfe_indexset;
	//!> constant above which the step width must always lie
	double constant_positivity;
	//!> constant to scale the linear interpolation
	double constant_interpolation;
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
	{
		std::vector<double> tempvector(x->size);
		minimizer->convertFromInternalType(x, tempvector);
		minimizer->functional.convertFromArrayType(
				tempvector,
				minimizer->value);
	}
	// evaluate the function (with array_type)
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
	{
		std::vector<double> tempvector(x->size);
		minimizer->convertFromInternalType(x, tempvector);
		minimizer->functional.convertFromArrayType(
				tempvector,
				minimizer->value);
	}
	// evaluate the gradient
	const T tmpgradient =
			minimizer->functional.gradient(minimizer->value);
	// convert gradient and store in g
	{
		std::vector<double> tempvector(g->size);
		minimizer->functional.convertToArrayType(
				tmpgradient,
				tempvector);
		minimizer->convertToInternalType(tempvector, g);
	}
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
	{
		std::vector<double> tempvector(x->size);
		minimizer->convertFromInternalType(x, tempvector);
		minimizer->functional.convertFromArrayType(
				tempvector,
				minimizer->value);
	}
	// evaluate the function and store in f
	*f = minimizer->functional(minimizer->value);
	// evaluate the gradient
	const T tmpgradient =
			minimizer->functional.gradient(minimizer->value);
	// convert gradient and store in g
	{
		std::vector<double> tempvector(g->size);
		minimizer->functional.convertToArrayType(
				tmpgradient,
				tempvector);
		minimizer->convertToInternalType(tempvector, g);
	}
}

//!> use this macro to make sure wrapper functions for gsl are instantiated
#define CONSTRUCT_FUNCTIONMINIMIZER(type) \
		template double FunctionMinimizer_FunctionCaller<type>(const gsl_vector *x,void *adata); \
		template void FunctionMinimizer_GradientCaller<type>(const gsl_vector *x,void *adata,gsl_vector *g); \
		template void FunctionMinimizer_FunctionGradientCaller<type>(const gsl_vector *x,void *adata,double *f,gsl_vector *g);

#include "Minimizations/Functions/FunctionMinimizer_impl.hpp"


#endif /* FUNCTIONMINIMIZER_HPP_ */
