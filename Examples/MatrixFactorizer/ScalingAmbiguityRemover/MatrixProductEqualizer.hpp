/*
 * MatrixProductEqualizer.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_SCALINGAMBIGUITYREMOVER_MATRIXPRODUCTEQUALIZER_HPP_
#define MATRIXFACTORIZERBASE_SCALINGAMBIGUITYREMOVER_MATRIXPRODUCTEQUALIZER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

/** Equalizes the scaling of both matrices by multiplying them with the
 * weighted average of both scaling factors over each individual scaling
 * factor.
 *
 * This is a sort of iterative way of obtaining rescaled matrices,
 * \sa MatrixProductRenormalizer
 */
struct MatrixProductEqualizer
{
	void operator()(
			Eigen::MatrixXd &_matrixone,
			Eigen::MatrixXd &_matrixtwo)
	{
		const double factorone = _matrixone.diagonal().maxCoeff();
		const double factortwo = _matrixtwo.diagonal().maxCoeff();
		double factor = .5*(factorone + factortwo);
		if (fabs(factor) > BASSOTOLERANCE) {
			if (factorone > factortwo)
				factor = 1./factor;
			_matrixone *= factor;
			_matrixtwo *= 1./factor;
		}
		//	if (_matrix.hasNaN())
		//		throw MinimizerIllegalNumber_exception()
		//		<< MinimizerIllegalNumber_variablename("matrix");
	}

};


#endif /* MATRIXFACTORIZERBASE_SCALINGAMBIGUITYREMOVER_MATRIXPRODUCTEQUALIZER_HPP_ */
