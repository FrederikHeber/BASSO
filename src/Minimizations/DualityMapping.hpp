/*
 * DualityMapping.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#ifndef DUALITYMAPPING_HPP_
#define DUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/LpNorm.hpp"

class DualityMappingUnitTest;

/** This class contains a duality mapping instance from a specific
 * lp space to its dual.
 *
 */
class DualityMapping
{
	//!> grant unit test DualityMappingUnitTest access to private parts.
	friend class DualityMappingUnitTest;
public:
	/** Constructor for class DualityMapping.
	 *
	 * \param _p p value of the used Lp norm
	 */
	DualityMapping(const double _p);
	~DualityMapping() {}

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
	Eigen::VectorXd operator()(
			const Eigen::VectorXd &_x,
			const double _power) const;

private:
	//!> value p of the Lp norm
	const double p;
	//!> tolerance for norm value
	mutable double tolerance;
	//!> lp norm object
	LpNorm lpnorm;
};

#endif /* DUALITYMAPPING_HPP_ */
