/*
 * GeneralMinimizer_DefaultDualityMappings.hpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#ifndef GENERALMINIMIZER_DEFAULTDUALITYMAPPINGS_HPP_
#define GENERALMINIMIZER_DEFAULTDUALITYMAPPINGS_HPP_

#include "BassoConfig.h"

#include "Minimizations/Mappings/LpDualityMapping.hpp"
#include "Minimizations/DualityMappings/DualityMappingsContainer.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

/** This equips the minimizer with the default differentiable duality mappings
 * and their inverse mappings.
 */
struct DefaultDualityMappings : public DualityMappingsContainer
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
		DualityMappingsContainer(
			_NormX,
			_PowerX,
			J_p,
			J_q),
		J_p(NormedSpaceFactory::DummySpace, val_NormX, _PowerX),
		J_q(NormedSpaceFactory::DummySpace, val_DualNormX, DualPowerX)
	{
		J_p.setTolerance(_tolerance);
		J_q.setTolerance(_tolerance);
	}

	//!> duality mapping object for space X
	const LpDualityMapping J_p;
	//!> duality mapping object for dual space X^*
	const LpDualityMapping J_q;
};


#endif /* GENERALMINIMIZER_DEFAULTDUALITYMAPPINGS_HPP_ */
