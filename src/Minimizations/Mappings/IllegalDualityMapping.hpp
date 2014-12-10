/*
 * IllegalDualityMapping.hpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#ifndef ILLEGALDUALITYMAPPING_HPP_
#define ILLEGALDUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Mappings/LpDualityMapping.hpp"

class IllegalDualityMappingUnitTest;

/** This class acts as a placeholder for a LpDualityMapping and ensures that this
 * function is never called in an algorithm.
 */
class IllegalDualityMapping :
	public LpDualityMapping
{
	//!> grant unit test access to private parts
	friend class IllegalDualityMappingUnitTest;
public:
	/** Constructor of class IllegalDualityMapping.
	 *
	 */
	IllegalDualityMapping();

	/** Getter for the source space of this mapping.
	 *
	 * @return ref to source space
	 */
	const NormedSpace_ptr_t& getSourceSpace() const;

	/** Getter for the target space of this mapping.
	 *
	 * @return ref to target space
	 */
	const NormedSpace_ptr_t& getTargetSpace() const;

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
	Mapping_ptr_t getAdjointMapping() const;

private:

	/** This function throws an assertion and must never be called.
	 * @param _x
	 * @param _power power of the duality mapping
	 * @return componentwise soft threshold of \a _x by \a _ lambda
	 */
	const Eigen::VectorXd operator()(
			const Eigen::VectorXd &_x) const;

};


#endif /* ILLEGALDUALITYMAPPING_HPP_ */
