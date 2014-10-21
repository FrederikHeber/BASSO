/*
 * GeneralMinimizer_RegularizedL1Norm.hpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#ifndef GENERALMINIMIZER_REGULARIZEDL1NORM_HPP_
#define GENERALMINIMIZER_REGULARIZEDL1NORM_HPP_

#include "BassoConfig.h"

#include "Minimizations/DualityMappings/DualityMappingsContainer.hpp"
#include "Minimizations/DualityMappings/IllegalDualityMapping.hpp"
#include "Minimizations/DualityMappings/SoftThresholdingOperator.hpp"

/** This equips the minimizer with a soft thresholding operator as required
 * by the regularized L1 norm ansatz.
 */
struct RegularizedL1Norm : public DualityMappingsContainer
{
	/** Constructor for class Regularized1Norm.
	 *
	 * @param _lambda value of the regularization parameter
	 */
	RegularizedL1Norm(
			const double _lambda) :
		DualityMappingsContainer(
			1.,	/* always have l1 norm here */
			1., /* this yields inf for q value */
			J_p,
			J_q),
			J_q(_lambda)
	{}

	//!> duality mapping object for space X
	const IllegalDualityMapping J_p;
	//!> duality mapping object for dual space X^*
	const SoftThresholdingOperator J_q;
};


#endif /* GENERALMINIMIZER_REGULARIZEDL1NORM_HPP_ */
