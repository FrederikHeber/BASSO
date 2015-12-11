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
};



#endif /* MINIMIZATIONS_FUNCTIONS_MINIMIZERS_FUNCTIONMINIMIZER_HPP_ */
