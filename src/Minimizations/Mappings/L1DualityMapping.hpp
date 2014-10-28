/*
 * L1DualityMapping.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef L1DUALITYMAPPING_HPP_
#define L1DUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/L1Norm.hpp"

class L1DualityMappingUnitTest;

/** This class contains a duality mapping instance from a specific
 * l_1 space to its dual.
 *
 */
class L1DualityMapping : public PowerTypeDualityMapping
{
	//!> grant unit test L1DualityMappingUnitTest access to private parts.
	friend class L1DualityMappingUnitTest;
public:
	/** Constructor for class L1DualityMapping.
	 *
	 */
	L1DualityMapping() {}

	/** Evaluates duality mapping at \a _x.
	 *
	 * \param _x point where to evaluate
	 * \param _power power of duality mapping's weight
	 */
	virtual const Eigen::VectorXd operator()(
			const Eigen::VectorXd &_x,
			const double _power) const;

protected:
	//!> lp norm object
	L1Norm l1norm;
};

#endif /* L1DUALITYMAPPING_HPP_ */
