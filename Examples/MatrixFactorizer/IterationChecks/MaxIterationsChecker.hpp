/*
 * MaxIterationsChecker.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_ITERATIONCHECKS_MAXITERATIONSCHECKER_HPP_
#define MATRIXFACTORIZERBASE_ITERATIONCHECKS_MAXITERATIONSCHECKER_HPP_

#include "BassoConfig.h"

/** Checks whether the maximum allowed iteration check as been exceeded.
 *
 */
struct MaxIterationsChecker
{
	/** Constructor for class MaxIterationsChecker
	 *
	 * @param _max_iterations maximum iterations allowed
	 */
	MaxIterationsChecker(
			const unsigned int _max_iterations
			) :
		max_iterations(_max_iterations)
	{}

	/** Functor that compares given \a _iterations against maximum
	 * allowed number
	 *
	 * @param _iterations current iterations done
	 * @return true - iteration count exceeded (stop)
	 */
	inline bool operator()(
			const unsigned int _iterations)
	{
		return _iterations > max_iterations;
	}

private:
	//!> internally stored maximum iterations count
	const unsigned int max_iterations;
};



#endif /* MATRIXFACTORIZERBASE_ITERATIONCHECKS_MAXITERATIONSCHECKER_HPP_ */
