/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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
