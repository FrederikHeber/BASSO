/*
 * LpDualityMapping.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LpDualityMapping.hpp"

#include <cmath>
#include <Eigen/Dense>

#include "Math/Helpers.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"

LpDualityMapping::LpDualityMapping(
		const double _p,
		const double _power) :
	PowerTypeDualityMapping(_power),
	p(_p),
	lpnorm(_p)
{}

LpDualityMapping::LpDualityMapping(
		const NormedSpace_ptr_t &_NormedSpaceRef,
		const double _p,
		const double _power) :
	PowerTypeDualityMapping(_NormedSpaceRef, _power),
	p(_p),
	lpnorm(_p)
{}

LpDualityMapping::LpDualityMapping(
		const double _p
		) :
	PowerTypeDualityMapping(_p),
	p(_p),
	lpnorm(_p)
{
	// we don't need to throw, norm will do this
}

LpDualityMapping::LpDualityMapping(
		const NormedSpace_ptr_t &_NormedSpaceRef,
		const double _p
		) :
	PowerTypeDualityMapping(_NormedSpaceRef, _p),
	p(_p),
	lpnorm(_p)
{
	// we don't need to throw, norm will do this
}

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
const Eigen::VectorXd LpDualityMapping::operator()(
		const Eigen::VectorXd &_x
		) const
{
	if (p == power) {
		// J=abs(x).^(p-1).*sign(x);
		Eigen::VectorXd Jx = _x.array().abs();
		for (int i=0;i<Jx.innerSize();++i)
			Jx[i] = ::pow(Jx[i], p - 1.) * Helpers::sign(_x[i]);
		return Jx;
	} else if (p < (double)power) {
		// J=norm(x,p)^(q-p)*abs(x).^(p-1).*sign(x);
		const double pnorm = ::pow(lpnorm(_x), (double)power-p);
		Eigen::VectorXd Jx = _x.array().abs();
		for (int i=0;i<Jx.innerSize();++i)
			Jx[i] = pnorm * ::pow(Jx[i], p - 1.) * Helpers::sign(_x[i]);
		return Jx;
	} else {
		const double norm = lpnorm(_x);
		if (norm < tolerance) {
			// J=zeros(size(x,1),1);
			Eigen::VectorXd Jx = Eigen::VectorXd::Zero(_x.innerSize());
			return Jx;
		} else {
			// J=n^(q-p)*abs(x).^(p-1).*sign(x);
			const double exponent = (double)power-p;
			const double pnorm = ::pow(norm, exponent);
			Eigen::VectorXd Jx = _x.array().abs();
			for (int i=0;i<Jx.innerSize();++i)
				Jx[i] = pnorm * ::pow(Jx[i], p - 1.) * Helpers::sign(_x[i]);
			return Jx;
		}
	}
}

Mapping_ptr_t LpDualityMapping::getAdjointMapping() const
{
	// calculate conjugate value
	const double q = Helpers::ConjugateValue(p);
	const double dualpower = Helpers::ConjugateValue(power);
	// and create instance with it
	PowerTypeDualityMapping_ptr_t instance(
			new LpDualityMapping(getTargetSpace(),q,dualpower));
	return instance;
}
