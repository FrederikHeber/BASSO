/*
 * L1DualityMapping.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "L1DualityMapping.hpp"

#include <cmath>

#include "Minimizations/Mappings/LInfinityDualityMapping.hpp"
#include "Minimizations/Norms/L1Norm.hpp"
#include "Math/Helpers.hpp"

/** General function to calculate the duality mapping.
 *
 *	In [Schöpfer et al., '06] some formulas for the duality mapping in
 *	Lp and other spaces are given. Note that also the power type of the
 *	duality mapping is denoted by p, making it ambigious.
 *
 *	We have settled on the following:
 *	The norm \f$ ||.|| \f$ is always the one of the space the argument
 *	of the duality mapping lives in, i.e. independent of the p in
 *	\f$ L_p \f$. We here refer to \f$ L_p \f$'s p as \a power, where
 *	\a p denotes the norm of the space.
 *
 *	With this convention the implementations in this function match
 *	the ones found in  [Schöpfer et al., '06].
 *
 * \param _x vector
 * \return dual element corresponding to one element of the duality mapping for
 * 		x
 */
const SpaceElement_ptr_t L1DualityMapping::operator()(
		const SpaceElement_ptr_t &_x
		) const
{
	// single-valued selection
	// J=norm(x,1)^(q-1)*sign(x);
	const Norm &l1norm = *getSourceSpace()->getNorm();

	const double factor = ::pow(l1norm(_x), (double)power-1.);
	SpaceElement_ptr_t sign_x = getTargetSpace()->createElement();
	*sign_x = _x->getSignVector()->getVectorRepresentation();
	*sign_x *= factor;
	return sign_x;
}

const Mapping_ptr_t L1DualityMapping::getAdjointMapping() const
{
	// calculate dual power
	const double dualpower = Helpers::ConjugateValue(power);
	// adjoint mapping is from l_infinity (from target to source, and
	// source is dual, where we assume spaces to be reflexive)
	PowerTypeDualityMapping_ptr_t instance(
			new LInfinityDualityMapping(
					getTargetSpace(), dualpower)
	);
	return instance;
}
