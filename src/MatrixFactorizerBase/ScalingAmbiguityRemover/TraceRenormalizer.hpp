/*
 * TraceRenormalizer.hpp
 *
 *  Created on: Jun 26, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZERBASE_SCALINGAMBIGUITYREMOVER_TRACERENORMALIZER_HPP_
#define MATRIXFACTORIZERBASE_SCALINGAMBIGUITYREMOVER_TRACERENORMALIZER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

/** Rescales the matrix by setting the maximum coefficient on its diagonal
 * to (plus oder minus) unity.
 */
struct TraceRenormalizer
{
	void operator()(Eigen::MatrixXd &_matrix)
	{
		const double factor = _matrix.diagonal().maxCoeff();
		if (fabs(factor) > BASSOTOLERANCE)
			_matrix *= 1./factor;
	//	if (_matrix.hasNaN())
	//		throw MinimizerIllegalNumber_exception()
	//		<< MinimizerIllegalNumber_variablename("matrix");
	}
};



#endif /* MATRIXFACTORIZERBASE_SCALINGAMBIGUITYREMOVER_TRACERENORMALIZER_HPP_ */
