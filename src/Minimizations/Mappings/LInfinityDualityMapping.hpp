/*
 * LInfinityDualityMapping.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef LINFINITYDUALITYMAPPING_HPP_
#define LINFINITYDUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/LInfinityNorm.hpp"

class LInfinityDualityMappingUnitTest;

/** This class contains a duality mapping instance from a specific
 * l_infty space to its dual.
 *
 */
class LInfinityDualityMapping : public PowerTypeDualityMapping
{
	//!> grant unit test LInfinityDualityMappingUnitTest access to private parts.
	friend class LInfinityDualityMappingUnitTest;
public:
	/** Constructor for class LInfinityDualityMapping.
	 *
	 */
	LInfinityDualityMapping() {}

	/** Evaluates duality mapping at \a _x.
	 *
	 * \param _x point where to evaluate
	 * \param _power power of duality mapping's weight
	 */
	virtual const Eigen::VectorXd operator()(
			const Eigen::VectorXd &_x,
			const double _power) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	PowerTypeDualityMapping_ptr_t getAdjointMapping() const;

protected:
	//!> lp norm object
	LInfinityNorm linftynorm;
};

#endif /* LINFINITYDUALITYMAPPING_HPP_ */
