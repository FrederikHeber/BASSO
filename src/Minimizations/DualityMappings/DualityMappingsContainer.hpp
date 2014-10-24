/*
 * DualityMappingsContainer.hpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#ifndef DUALITYMAPPINGSCONTAINER_HPP_
#define DUALITYMAPPINGSCONTAINER_HPP_

#include "BassoConfig.h"

#include "Minimizations/Mappings/LpDualityMapping.hpp"

/** This class defines a container for the duality mappings for a single
 * space X, to be used for example by the GeneralMinimizer class.
 */
struct DualityMappingsContainer
{
	DualityMappingsContainer(
			const double _NormX,
			const double _PowerX,
			const LpDualityMapping &_J_p,
			const LpDualityMapping &_J_q) :
		val_NormX(_NormX),
		val_DualNormX(val_NormX/(val_NormX - 1.)),
		PowerX(_PowerX),
		DualPowerX(PowerX/(PowerX - 1.)),
		J_p(_J_p),
		J_q(_J_q)
	{}

	DualityMappingsContainer(
			const DualityMappingsContainer &_container) :
		val_NormX(_container.val_NormX),
		val_DualNormX(_container.val_DualNormX),
		PowerX(_container.PowerX),
		DualPowerX(_container.DualPowerX),
		J_p(_container.J_p),
		J_q(_container.J_q)
	{}

	//!> Lp norm of space X: p
	const double val_NormX;
	//!> Lp norm of dual space to X: q
	const double val_DualNormX;
	//!> power of dual map J_p
	const double PowerX;
	//!> power of dual map J_q
	const double DualPowerX;

	//!> reference to duality mapping object for space X
	const LpDualityMapping &J_p;
	//!> reference to duality mapping object for dual space X^*
	const LpDualityMapping &J_q;
};


#endif /* DUALITYMAPPINGSCONTAINER_HPP_ */
