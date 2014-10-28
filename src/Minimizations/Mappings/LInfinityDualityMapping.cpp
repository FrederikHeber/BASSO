/*
 * LInfinityDualityMapping.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LInfinityDualityMapping.hpp"

#include <cmath>
#include <Eigen/Dense>

#include "Minimizations/Mappings/L1DualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/MinimizationExceptions.hpp"
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
const Eigen::VectorXd LInfinityDualityMapping::operator()(
		const Eigen::VectorXd &_x
		) const
{
	// [xNorm,k]=max(abs(x));
	unsigned int rowMax;
	unsigned int colMax;
	double factor = _x.array().abs().maxCoeff(&rowMax, &colMax);
	factor = ::pow(factor, (double)power-1.) * Helpers::sign(_x[rowMax]);
	// J=xNorm^(q-1)*sign(x(k,1))*circshift(eye(size(x)),[k-1 0]);
	Eigen::VectorXd temp = Eigen::VectorXd::Zero(_x.innerSize());
	temp[0] = 1.;
	const Eigen::VectorXd Jx = Helpers::circshift(temp, rowMax); // no -1 here, as index starts at 0 here, not 1
	return Jx * factor;
}

Mapping_ptr_t LInfinityDualityMapping::getAdjointMapping() const
{
	// calculate dual power
	const double dualpower = Helpers::ConjugateValue(power);
	// adjoint mapping is from l_1
	PowerTypeDualityMapping_ptr_t instance(
			new L1DualityMapping(dualpower)
	);
	return instance;
}
