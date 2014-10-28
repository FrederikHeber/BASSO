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

/** This class acts as a placeholder for a LpDualityMapping and ensures that this
 * function is never called in an algorithm.
 */
class IllegalDualityMapping :
	public LpDualityMapping
{
public:
	/** Constructor of class IllegalDualityMapping.
	 *
	 */
	IllegalDualityMapping();

	/** This function throws an assertion and must never be called.
	 * @param _x
	 * @param _power power of the duality mapping
	 * @return componentwise soft threshold of \a _x by \a _ lambda
	 */
	const Eigen::VectorXd operator()(
			const Eigen::VectorXd &_x) const;
};


#endif /* ILLEGALDUALITYMAPPING_HPP_ */
