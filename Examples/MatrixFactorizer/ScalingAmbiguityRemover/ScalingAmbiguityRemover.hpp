/*
 * ScalingAmbiguityRemover.hpp
 *
 *  Created on: Nov 23, 2015
 *      Author: heber
 */

#ifndef SOLVERS_SCALINGAMBIGUITYREMOVER_HPP_
#define SOLVERS_SCALINGAMBIGUITYREMOVER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

/** This class defines the interface for a scaling ambiguity
 * remover.
 */
class ScalingAmbiguityRemover
{
	virtual double operator()(
			Eigen::MatrixXd &_matrixone,
			Eigen::MatrixXd &_matrixtwo) const = 0;
};


#endif /* SOLVERS_SCALINGAMBIGUITYREMOVER_HPP_ */
