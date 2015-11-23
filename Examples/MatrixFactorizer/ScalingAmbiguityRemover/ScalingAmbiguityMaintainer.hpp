/*
 * ScalingAmbiguityMaintainer.hpp
 *
 *  Created on: Nov 23, 2015
 *      Author: heber
 */

#ifndef SOLVERS_SCALINGAMBIGUITYMAINTAINER_HPP_
#define SOLVERS_SCALINGAMBIGUITYMAINTAINER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "ScalingAmbiguityRemover.hpp"

/** This implementation of a ScalingAmbiguityRemover does nothing at
 * all.
 */
struct ScalingAmbiguityMaintainer : public ScalingAmbiguityRemover
{
	double operator()(
			Eigen::MatrixXd &_matrixone,
			Eigen::MatrixXd &_matrixtwo) const
	{
		const double factorone = _matrixone.diagonal().maxCoeff();
		const double factortwo = _matrixtwo.diagonal().maxCoeff();
		double factor = .5*(factorone + factortwo);
		return factor;
	}
};


#endif /* SOLVERS_SCALINGAMBIGUITYMAINTAINER_HPP_ */
