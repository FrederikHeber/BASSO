/*
 * SoftThresholdingOperator.hpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#ifndef SOFTTHRESHOLDINGOPERATOR_HPP_
#define SOFTTHRESHOLDINGOPERATOR_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/DualityMappings/DualityMapping.hpp"

class SoftThresholdingOperator :
	public DualityMapping
{
public:
	/** Default constructor of class SoftThresholdingOperator.
	 *
	 */
	SoftThresholdingOperator();

	/** Default constructor of class SoftThresholdingOperator.
	 *
	 * @param _lambda soft thresholding parameter
	 */
	SoftThresholdingOperator(const double _lambda);

	/** Evaluates for the given \a _x the soft thresholding result with respect
	 * to \a _lambda.
	 * @param _x vector to soft-threshold
	 * @param _power \b not used, kept for signature compatibility
	 * @return componentwise soft threshold of \a _x by \a lambda
	 */
	Eigen::VectorXd operator()(
			const Eigen::VectorXd &_x,
			const double) const;

	/** Setter for soft thresholding parameter \a lambda.
	 *
	 * @param _lambda new value for parameter
	 */
	void setLambda(const double _lambda) const
	{
		const_cast<double &>(lambda) = _lambda;
	}

private:
	//!> soft thresholding parameter
	const double lambda;
};


#endif /* SOFTTHRESHOLDINGOPERATOR_HPP_ */