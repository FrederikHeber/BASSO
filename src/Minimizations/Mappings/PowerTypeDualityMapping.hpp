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

#include "Minimizations/types.hpp"
#include "Minimizations/Mappings/DualityMapping.hpp"

/** This class defines mappings from a space to its dual space with
 * elements being the gradient of the norm of the source space.
 */
class PowerTypeDualityMapping : public DualityMapping
{
public:
	/** Default constructor for class PowerTypeDualityMapping.
	 *
	 * @param _power power type of this duality mapping
	 */
	PowerTypeDualityMapping(
			const double _power) :
		power(_power),
		tolerance(BASSOTOLERANCE)
	{}

	/** Default constructor for class PowerTypeDualityMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _power power type of this duality mapping
	 */
	PowerTypeDualityMapping(
			const NormedSpace_ptr_t &_NormedSpaceRef,
			const double _power) :
		DualityMapping(_NormedSpaceRef),
		power(_power),
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
			const Eigen::VectorXd &_x) const = 0;

	/** Mapping function.
	 *
	 * @param _sourceelement element to map/transform
	 * @return new transformed/mapped element
	 */
	SpaceElement_ptr_t operator()(
			const SpaceElement_ptr_t &_sourceelement
			) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	virtual Mapping_ptr_t getAdjointMapping() const = 0;

	/** Returns the tolerance.
	 *
	 * @return tolerance tolerance threshold
	 */
	const double getTolerance() const
	{ return tolerance; }

	/** Returns the power type of th mapping's weight function.
	 *
	 * @return power type
	 */
	const double getPower() const
	{ return power; }

protected:
	//!> power type of the weight function of this duality mapping
	const double power;
	//!> tolerance for norm value
	mutable double tolerance;
};


#endif /* POWERTYPEDUALITYMAPPING_HPP_ */
