/*
 * LpDualityMapping.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#ifndef LPDUALITYMAPPING_HPP_
#define LPDUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Norms/LpNorm.hpp"

class LpDualityMappingUnitTest;

/** This class contains a duality mapping instance from a specific
 * lp space to its dual.
 *
 */
class LpDualityMapping
{
	//!> grant unit test LpDualityMappingUnitTest access to private parts.
	friend class LpDualityMappingUnitTest;
public:
	/** Constructor for class LpDualityMapping.
	 *
	 * \param _p p value of the used Lp norm
	 */
	LpDualityMapping(const double _p);
	virtual ~LpDualityMapping() {}

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
			const double _power) const;

protected:
	//!> value p of the Lp norm
	const double p;
	//!> tolerance for norm value
	mutable double tolerance;
	//!> lp norm object
	LpNorm lpnorm;
};

#endif /* LPDUALITYMAPPING_HPP_ */
