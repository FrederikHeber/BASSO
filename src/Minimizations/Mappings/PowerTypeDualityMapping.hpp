/*
 * PowerTypeDualityMapping.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef POWERTYPEDUALITYMAPPING_HPP_
#define POWERTYPEDUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

/** This class defines mappings from a space to its dual space with
 * elements being the gradient of the norm of the source space.
 */
class PowerTypeDualityMapping
{
public:
	/** Default constructor for class PowerTypeDualityMapping.
	 *
	 */
	PowerTypeDualityMapping() :
		tolerance(BASSOTOLERANCE)
	{}

	/** Setter for internal tolerance.
	 *
	 * \param _tolerance value to set to
	 */
	void setTolerance(const double _tolerance) const
	{ tolerance = _tolerance; }

	/** Evaluates duality mapping at \a _x.
	 *
	 * \param _x point where to evaluate
	 * \param _power power of duality mapping's weight
	 */
	virtual const Eigen::VectorXd operator()(
			const Eigen::VectorXd &_x,
			const double _power) const = 0;

	/** Returns the tolerance.
	 *
	 * @return tolerance tolerance threshold
	 */
	const double getTolerance() const
	{ return tolerance; }

protected:
	//!> tolerance for norm value
	mutable double tolerance;
};


#endif /* POWERTYPEDUALITYMAPPING_HPP_ */
