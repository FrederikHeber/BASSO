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
		const NormedSpace_ptr_t &_NormedSpaceRef,
		const double _power) :
	PowerTypeDualityMapping(_NormedSpaceRef, _power)
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
const SpaceElement_ptr_t LpDualityMapping::operator()(
		const SpaceElement_ptr_t &_x
		) const
{
	const Norm &lpnorm = *getSourceSpace()->getNorm();
	const double p = lpnorm.getPvalue();
	if (p == power) {
		// J=abs(x).^(p-1).*sign(x);
		SpaceElement_ptr_t Jx =
				ElementCreator::create(
						getTargetSpace(),
						RepresentationAdvocate::get(_x->getAbsVector()));
		SpaceElement_ptr_t sign_x = _x->getSignVector();
		for (unsigned int i=0;i<Jx->getSpace()->getDimension();++i)
			(*Jx)[i] = ::pow((*Jx)[i], p - 1.) * (*sign_x)[i];
		return Jx;
	} else if (p < (double)power) {
		// J=norm(x,p)^(q-p)*abs(x).^(p-1).*sign(x);
		const double pnorm = ::pow(lpnorm(_x), (double)power-p);
		SpaceElement_ptr_t Jx =
				ElementCreator::create(
						getTargetSpace(),
						RepresentationAdvocate::get(_x->getAbsVector()));
		SpaceElement_ptr_t sign_x = _x->getSignVector();
		for (unsigned int i=0;i<Jx->getSpace()->getDimension();++i)
			(*Jx)[i] = pnorm * ::pow((*Jx)[i], p - 1.) * (*sign_x)[i];
		return Jx;
	} else {
		const double norm = lpnorm(_x);
		if (norm < tolerance) {
			// J=zeros(size(x,1),1);
			SpaceElement_ptr_t Jx = getTargetSpace()->createElement();
			return Jx;
		} else {
			// J=n^(q-p)*abs(x).^(p-1).*sign(x);
			const double exponent = (double)power-p;
			const double pnorm = ::pow(norm, exponent);
			SpaceElement_ptr_t Jx =
					ElementCreator::create(
							getTargetSpace(),
							RepresentationAdvocate::get(_x->getAbsVector()));
			SpaceElement_ptr_t sign_x = _x->getSignVector();
			for (unsigned int i=0;i<Jx->getSpace()->getDimension();++i)
				(*Jx)[i] = pnorm * ::pow((*Jx)[i], p - 1.) * (*sign_x)[i];
			return Jx;
		}
	}
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
