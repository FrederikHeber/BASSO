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
	 * @param _val_NormX p value of the Lp space
	 * @param _PowerX power of the Lp duality mapping
	 * @param _tolerance tolerance value of the duality mappings
	 */
	RegularizedL1Norm(
			const double _NormX,
			const double _PowerX) :
		DualityMappingsContainer(
			_NormX,
			_PowerX,
			J_p,
			J_q)
	{}

	//!> duality mapping object for space X
	const IllegalDualityMapping J_p;
	//!> duality mapping object for dual space X^*
	const SoftThresholdingOperator J_q;
};


#endif /* GENERALMINIMIZER_REGULARIZEDL1NORM_HPP_ */
