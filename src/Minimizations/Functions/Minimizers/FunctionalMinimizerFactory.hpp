/*
 * FunctionalMinimizerFactory.hpp
 *
 *  Created on: Dec 11, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_FUNCTIONS_MINIMIZERS_FUNCTIONALMINIMIZERFACTORY_HPP_
#define MINIMIZATIONS_FUNCTIONS_MINIMIZERS_FUNCTIONALMINIMIZERFACTORY_HPP_

#include "BassoConfig.h"

#include <string>

#include "FunctionalMinimizer.hpp"

/** Factory pattern for the Minimizer of a function.
 *
 * This pattern is property-based, i.e. small cstor, lots of setters.
 *
 */
struct FunctionalMinimizerFactory
{
	//!> enumeration of all available minimization libraries (for line search)
	enum MinimizationLibraries {
		gnuscientificlibrary=0,
		nonlinearoptimization=1,
		nonconvex_regularizedl1=2,
		MAX_MinimizationLibraries
	};

	/** Creates a minimizer instance for the current parameters.
	 *
	 * @param _N number of search directions
	 * @param _functional the function to minimize
	 * @param _value the external value as workspace for minimization
	 * @return created instance
	 */
	template <class T>
	static typename FunctionalMinimizer<T>::ptr_t create(
			const unsigned int _N,
			const MinimizationFunctional<T> &_functional,
			T &_value
			);

	/** Creates a minimizer instance for the current parameters, with inexact
	 * line search according to Wolfe conditions
	 *
	 * @param _N number of search directions
	 * @param _functional the function to minimize
	 * @param _value the external value as workspace for minimization
	 * @param _constant_positivity positivity constant for Wolfe conditions
	 * @param _Wolfe_indexset search directions on which to apply search directions
	 * @param _constant_interpolation Wolfe constant for stronger than linear
	 * @return created instance
	 */
	template <class T>
	static typename FunctionalMinimizer<T>::ptr_t create(
			const unsigned int _N,
			const MinimizationFunctional<T> &_functional,
			T &_value,
			const double _constant_positivity,
			const std::vector<unsigned int> &_Wolfe_indexset,
			const double _constant_interpolation
			);

	/** Setter for MinLib via \a _name as string.
	 *
	 * @param _name name of library to use for minimization
	 */
	static void setMinLib(const std::string &_name);

	/** Checks whether \a _name represents a valid name for a
	 * minimization library.
	 *
	 * @param _name name of library to check
	 * @return true - valid name, false - library name unknown
	 */
	static bool isValidName(const std::string &_name);

private:
	//!> name of each known instance type (Don't forget to add enum to InstanceType)
	static const std::string TypeNames[];

	//!> produced type of minimization instance
	static enum MinimizationLibraries CurrentMinLib;

	//!> maximum number of iterations
	static unsigned int maxiterations;
};

#include "FunctionalMinimizerFactory_impl.hpp"


#endif /* MINIMIZATIONS_FUNCTIONS_MINIMIZERS_FUNCTIONALMINIMIZERFACTORY_HPP_ */
