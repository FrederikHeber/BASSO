/*
 * RelativeChangeResidualChecker.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_ITERATIONCHECKS_RELATIVECHANGERESIDUALCHECKER_HPP_
#define MATRIXFACTORIZERBASE_ITERATIONCHECKS_RELATIVECHANGERESIDUALCHECKER_HPP_

#include "BassoConfig.h"

/** Checks the relative residual against a stored threshold.
 *
 */
struct RelativeChangeResidualChecker
{
	/** Constructor for class RelativeResidualChecker.
	 *
	 * @param _residual residual threshold to compare against.
	 */
	RelativeChangeResidualChecker(
			const double _delta
			) :
		delta(_delta),
		oldresidual(0.)
	{}

	/** Functor for checking the residual.
	 *
	 * @param _residual residual to check
	 * @return true - stop, relative change in residual is below threshold
	 */
	inline bool operator()(
			const double _residual)
	{
		bool result;
		if (fabs(_residual) > BASSOTOLERANCE)
			result = (fabs(oldresidual -_residual)/_residual) < delta;
		else
			result = true;
		oldresidual = _residual;
		return result;
	}

private:
	//!> internally stored residual threshold
	const double delta;

	//!> contains residual from last iteration
	double oldresidual;
};


#endif /* MATRIXFACTORIZERBASE_ITERATIONCHECKS_RELATIVECHANGERESIDUALCHECKER_HPP_ */
