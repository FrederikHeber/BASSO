/*
 * L1DualityMapping.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "L1DualityMapping.hpp"

#include <cmath>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Mappings/Specifics/LInfinityDualityMapping.hpp"
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
 * \param _Jx duality mapped \a _x
 */
void L1DualityMapping::operator()(
		const SpaceElement_ptr_t &_x,
		SpaceElement_ptr_t &_Jx
		) const
{
	// start timing
	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	// single-valued selection
	// J=norm(x,1)^(q-1)*sign(x);
	assert( getSourceSpace().get() == _x->getSpace().get() );
	assert( getTargetSpace().get() == _Jx->getSpace().get() );
	const Norm &l1norm = *_x->getSpace()->getNorm();

	const double factor = ::pow(l1norm(_x), (double)power-1.);
	assert( _Jx->getSpace()->getDimension() == _x->getSpace()->getDimension() );
	const Eigen::VectorXd &vector = RepresentationAdvocate::get(_x);
	RepresentationAdvocate::set(_Jx,
			factor * vector.array().cwiseProduct(
					vector.array().abs().cwiseInverse()));

	// finish timing
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	timing += timing_end - timing_start;
	++count;
}

void L1DualityMapping::getMinimumInfimum(
		const SpaceElement_ptr_t &_x,
		const SpaceElement_ptr_t &_y,
		SpaceElement_ptr_t &_Jx) const
{
	// l1 is (component-wise) only not differentiable where we intersect with axis
	// hence, for all intersections check whether flipped sign is better
	assert( getSourceSpace().get() == _x->getSpace().get() );
	assert( getTargetSpace().get() == _Jx->getSpace().get() );
	operator()(_x, _Jx);
	for (unsigned int i=0;i<_x->getSpace()->getDimension();++i) {
		if (fabs((*_x)[i]) < BASSOTOLERANCE ) {
			if ((*_Jx)[i] * (*_y)[i] > 0)
				(*_Jx)[i] *= -1.;
		}
	}
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
