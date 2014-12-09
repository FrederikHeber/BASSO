/*
 * LInfinityDualityMapping.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LInfinityDualityMapping.hpp"

#include <cmath>

#include "Minimizations/Mappings/Specifics/L1DualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Math/Helpers.hpp"

/** General function to calculate the duality mapping.
 *
 *	In [Schöpfer et al., '06] some formulas for the duality mapping in
 *	Lp and other spaces are given. Note that also the power type of the
 *	duality mapping is denoted by p, making it ambiguous.
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
const SpaceElement_ptr_t LInfinityDualityMapping::operator()(
		const SpaceElement_ptr_t &_x
		) const
{
	// start timing
	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	// [xNorm,k]=max(abs(x));
	const std::pair<double,int> factor_index =
			_x->getMaxCoefficientAndIndex();
	const double factor = ::pow(factor_index.first, (double)power-1.)
			* Helpers::sign((*_x)[factor_index.second]);
	// J=xNorm^(q-1)*sign(x(k,1))*circshift(eye(size(x)),[k-1 0]);
	SpaceElement_ptr_t temp = getTargetSpace()->createElement();
	(*temp)[0] = 1.;
	SpaceElement_ptr_t Jx =
			temp->getCircShiftedVector(factor_index.second);  // no -1 here, as index starts at 0 here, not 1
	*Jx *= factor;

	// finish timing
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	timing += timing_end - timing_start;
	++count;

	return Jx;
}

const Mapping_ptr_t LInfinityDualityMapping::getAdjointMapping() const
{
	// calculate dual power
	const double dualpower = Helpers::ConjugateValue(power);
	// adjoint mapping is from l_1
	PowerTypeDualityMapping_ptr_t instance(
			new L1DualityMapping(getTargetSpace(),dualpower)
	);
	return instance;
}
