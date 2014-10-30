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

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"

class LpDualityMappingUnitTest;

/** This class contains a duality mapping instance from a specific
 * l_p space with 1 < p < infty to its dual.
 *
 */
class LpDualityMapping : public PowerTypeDualityMapping
{
	//!> grant unit test LpDualityMappingUnitTest access to private parts.
	friend class LpDualityMappingUnitTest;
public:
	/** Constructor for class LpDualityMapping.
	 *
	 * \param _p p value of the used Lp norm
	 * @param _power power type of this duality mapping
	 */
	LpDualityMapping(
			const double _p,
			const double _power);

	/** Constructor for class LpDualityMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * \param _p p value of the used Lp norm
	 * @param _power power type of this duality mapping
	 */
	LpDualityMapping(
			const NormedSpace_ptr_t &_NormedSpaceRef,
			const double _p,
			const double _power);

	/** Constructor for class LpDualityMapping.
	 *
	 * Here, we use \a _p also for the power of the weight function.
	 *
	 * \param _p p value of the used Lp norm
	 */
	LpDualityMapping(const double _p);

	/** Constructor for class LpDualityMapping.
	 *
	 * Here, we use \a _p also for the power of the weight function.
	 *
	 * @param _NormedSpaceRef reference to space
	 * \param _p p value of the used Lp norm
	 */
	LpDualityMapping(
			const NormedSpace_ptr_t &_NormedSpaceRef,
			const double _p);

	/** Evaluates duality mapping at \a _x.
	 *
	 * \param _x point where to evaluate
	 */
	virtual const Eigen::VectorXd operator()(
			const Eigen::VectorXd &_x) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	Mapping_ptr_t getAdjointMapping() const;

protected:
	//!> value p of the Lp norm
	const double p;
	//!> lp norm object
	LpNorm lpnorm;
};

#endif /* LPDUALITYMAPPING_HPP_ */
