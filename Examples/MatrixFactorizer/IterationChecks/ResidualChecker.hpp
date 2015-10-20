/*
 * ResidualChecker.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_ITERATIONCHECKS_RESIDUALCHECKER_HPP_
#define MATRIXFACTORIZERBASE_ITERATIONCHECKS_RESIDUALCHECKER_HPP_

#include "BassoConfig.h"

/** Checks whether the residual has dropped below some threshold.
 *
 */
struct ResidualChecker
{
	/** Constructor for class ResidualChecker.
	 *
	 * @param _delta delta threshold to compare against.
	 */
	ResidualChecker(
			const double _delta
			) :
		delta(_delta)
	{}

	/** Functor for checking the residual.
	 *
	 * @param _residual current residual
	 * @return true - residual is less than given threshold \a delta
	 */
	inline bool operator()(
			const double _residual)
	{
		return _residual < delta;
	}

private:
	//!> internally stored residual threshold
	const double delta;
};


#endif /* MATRIXFACTORIZERBASE_ITERATIONCHECKS_RESIDUALCHECKER_HPP_ */
