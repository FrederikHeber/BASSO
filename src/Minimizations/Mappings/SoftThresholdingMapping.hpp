/*
 * SoftThresholdingMapping.hpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#ifndef SOFTTHRESHOLDINGMAPPING_HPP_
#define SOFTTHRESHOLDINGMAPPING_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Mappings/L1DualityMapping.hpp"

class SoftThresholdingMapping :
	public L1DualityMapping
{
public:
	/** Default constructor of class SoftThresholdingMapping.
	 *
	 */
	SoftThresholdingMapping();

	/** Default constructor of class SoftThresholdingMapping.
	 *
	 * @param _lambda soft thresholding parameter
	 */
	SoftThresholdingMapping(const double _lambda);

	/** Evaluates for the given \a _x the soft thresholding result with respect
	 * to \a _lambda.
	 * @param _x vector to soft-threshold
	 * @return componentwise soft threshold of \a _x by \a lambda
	 */
	const Eigen::VectorXd operator()(
			const Eigen::VectorXd &_x) const;

	/** Setter for soft thresholding parameter \a lambda.
	 *
	 * @param _lambda new value for parameter
	 */
	void setLambda(const double _lambda) const
	{
		const_cast<double &>(lambda) = _lambda;
	}

	/** Getter for the regularization parameter.
	 *
	 * @return regularization parameter
	 */
	const double getLambda() const
	{ return lambda; }

private:
	//!> soft thresholding parameter
	const double lambda;
};


#endif /* SOFTTHRESHOLDINGMAPPING_HPP_ */
