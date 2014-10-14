/*
 * GeneralMinimizer_RegularizedL1Norm.hpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#ifndef GENERALMINIMIZER_REGULARIZEDL1NORM_HPP_
#define GENERALMINIMIZER_REGULARIZEDL1NORM_HPP_

#include "Minimizations/DualityMappings/IllegalDualityMapping.hpp"
#include "Minimizations/DualityMappings/SoftThresholdingOperator.hpp"

/** This equips the minimizer with a soft thresholding operator as required
 * by the regularized L1 norm ansatz.
 */
struct RegularizedL1Norm
{
	/** Constructor for class Regularized1Norm.
	 *
	 * @param _val_NormX p value of the Lp space
	 * @param _PowerX power of the Lp duality mapping
	 * @param _tolerance tolerance value of the duality mappings
	 */
	RegularizedL1Norm(
			const double _NormX,
			const double _PowerX,
			const double _tolerance) :
		val_NormX(_NormX),
		val_DualNormX(val_NormX/(val_NormX - 1.)),
		PowerX(_PowerX),
		DualPowerX(PowerX/(PowerX - 1.))
	{}

	//!> Lp norm of space X: p
	const double val_NormX;
	//!> Lp norm of dual space to X: q
	const double val_DualNormX;
	//!> power of dual map J_p
	const double PowerX;
	//!> power of dual map J_q
	const double DualPowerX;

	//!> duality mapping object for space X
	const IllegalDualityMapping J_p;
	//!> duality mapping object for dual space X^*
	const SoftThresholdingOperator J_q;
};


#endif /* GENERALMINIMIZER_REGULARIZEDL1NORM_HPP_ */
