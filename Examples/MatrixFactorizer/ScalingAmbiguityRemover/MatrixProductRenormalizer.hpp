/*
 * MatrixProductRenormalizer.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_SCALINGAMBIGUITYREMOVER_MATRIXPRODUCTRENORMALIZER_HPP_
#define MATRIXFACTORIZERBASE_SCALINGAMBIGUITYREMOVER_MATRIXPRODUCTRENORMALIZER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "ScalingAmbiguityRemover.hpp"

/** Calculates the square root of the product of either matrix factors and
 * rescales by it over the original factor.
 */
struct MatrixProductRenormalizer : public ScalingAmbiguityRemover
{
	double operator()(
			Eigen::MatrixXd &_matrixone,
			Eigen::MatrixXd &_matrixtwo) const
	{
		double maxfactor = 0.;
		for (unsigned int dim = 0; dim <_matrixone.outerSize();++dim ) {
			const double factor = _matrixone.col(dim).lpNorm<Eigen::Infinity>();
			if (fabs(factor) > BASSOTOLERANCE) {
				_matrixone.col(dim) *= 1./factor;
				_matrixtwo.row(dim) *= factor;
			}
			maxfactor = std::max(factor, maxfactor);
		}

		return maxfactor;
	//	if (_matrix.hasNaN())
	//		throw MinimizerIllegalNumber_exception()
	//		<< MinimizerIllegalNumber_variablename("matrix");
	}

};

#endif /* MATRIXFACTORIZERBASE_SCALINGAMBIGUITYREMOVER_MATRIXPRODUCTRENORMALIZER_HPP_ */
