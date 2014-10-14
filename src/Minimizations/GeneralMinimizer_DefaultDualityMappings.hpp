/*
 * GeneralMinimizer_DefaultDualityMappings.hpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#ifndef GENERALMINIMIZER_DEFAULTDUALITYMAPPINGS_HPP_
#define GENERALMINIMIZER_DEFAULTDUALITYMAPPINGS_HPP_

#include "Minimizations/DualityMappings/DualityMapping.hpp"

/** This equips the minimizer with the default differentiable duality mappings
 * and their inverse mappings.
 */
struct DefaultDualityMappings
{
	/** Constructor for class DefaultDualityMappings.
	 *
	 * q value is calculated as conjugated to given \a _NormX and \a _PowerX.
	 *
	 * @param _val_NormX p value of the Lp space
	 * @param _PowerX power of the Lp duality mapping
	 * @param _tolerance tolerance value of the duality mappings
	 */
	DefaultDualityMappings(
			const double _NormX,
			const double _PowerX,
			const double _tolerance) :
		val_NormX(_NormX),
		val_DualNormX(val_NormX/(val_NormX - 1.)),
		PowerX(_PowerX),
		DualPowerX(PowerX/(PowerX - 1.)),
		J_p(val_NormX),
		J_q(val_DualNormX)
	{
		J_p.setTolerance(_tolerance);
		J_q.setTolerance(_tolerance);
	}

	//!> Lp norm of space X: p
	const double val_NormX;
	//!> Lp norm of dual space to X: q
	const double val_DualNormX;
	//!> power of dual map J_p
	const double PowerX;
	//!> power of dual map J_q
	const double DualPowerX;

	//!> duality mapping object for space X
	const DualityMapping J_p;
	//!> duality mapping object for dual space X^*
	const DualityMapping J_q;
};


#endif /* GENERALMINIMIZER_DEFAULTDUALITYMAPPINGS_HPP_ */
