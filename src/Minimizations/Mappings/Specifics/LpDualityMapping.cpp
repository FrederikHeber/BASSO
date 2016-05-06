/*
 * LpDualityMapping.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LpDualityMapping.hpp"

#include <cmath>

#include "Math/Helpers.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"

LpDualityMapping::LpDualityMapping(
		const NormedSpace_weakptr_t &_NormedSpaceRef,
		const double _power) :
	PowerTypeDualityMapping(_NormedSpaceRef, _power),
	count(0),
	timing(boost::chrono::nanoseconds(0))
{}

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
void LpDualityMapping::operator()(
		const SpaceElement_ptr_t &_x,
		SpaceElement_ptr_t &_Jx
		) const
{
	// start timing
	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	const Norm &lpnorm = *getSourceSpace()->getNorm();
	const double p = lpnorm.getPvalue();
	// (also in case: p == power)
	// J=abs(x).^(p-1).*sign(x);
	// need to do it this complicatedly as spaces aren't right
	RepresentationAdvocate::set(
			_Jx,
			RepresentationAdvocate::get(_x).array().abs());
	_Jx->pow(p-1.);
	RepresentationAdvocate::set(
			_Jx,
			RepresentationAdvocate::get(_Jx).cwiseProduct(
					RepresentationAdvocate::get(_x).unaryExpr(std::ptr_fun(Helpers::sign)))
	);
	if (p != power) {
		const double norm = lpnorm(_x);
		const double pnorm = ::pow(norm, (double)power-p);
		if ((p < (double)power) || (norm > tolerance)) {
			// J=norm(x,p)^(q-p)*abs(x).^(p-1).*sign(x);
			*_Jx *= pnorm;
		} else {
			// J=zeros(size(x,1),1);
			_Jx->setZero();
		}
	}

	// finish timing
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	timing += timing_end - timing_start;
	++count;
}

const Mapping_ptr_t LpDualityMapping::getAdjointMapping() const
{
	// calculate conjugate value
	const double dualpower = Helpers::ConjugateValue(power);
	// and create instance with it
	PowerTypeDualityMapping_ptr_t instance(
			new LpDualityMapping(getTargetSpace(),dualpower));
	return instance;
}
