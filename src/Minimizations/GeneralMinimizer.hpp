/*
 * GeneralMinimizer.hpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#ifndef GENERALMINIMIZER_HPP_
#define GENERALMINIMIZER_HPP_

#include <Eigen/Dense>

/** This class describes the interface to a general minimizer.
 *
 */
class GeneralMinimizer
{
public:
	GeneralMinimizer() {}
	virtual ~GeneralMinimizer() {}

	/** Internal structure for return values.
	 *
	 */
	struct ReturnValues
	{
		//!> solution vector
		Eigen::VectorXd solution;
		//!> remaining residuum
		double residuum;
		//!> number of outer iterations till solution
		int NumberOuterIterations;
	};


};


#endif /* GENERALMINIMIZER_HPP_ */
