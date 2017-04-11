/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
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
 * BregmanProjectionFunctional.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "BregmanProjectionFunctional.hpp"

#include <cmath>
#include <Eigen/Dense>
#include <numeric>

#include "Log/Logging.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Norms/Norm.hpp"

BregmanProjectionFunctional::BregmanProjectionFunctional(
		const Norm &_dualnorm,
		const Mapping &_J_q,
		const double _dualpower,
		const std::vector<SpaceElement_ptr_t> &_U,
		const std::vector<double> &_alpha
		) :
	dualpower(_dualpower),
	dualnorm(_dualnorm),
	J_q(_J_q),
	U(_U),
	alpha(_alpha),
	resx((*U.begin())->getSpace()->createElement()),
	dual_resx(resx->getSpace()->getDualSpace()->createElement()),
	zeroVec((*U.begin())->getSpace()->createElement())
{
	assert ( U.size() == alpha.size() );
	assert( !U.empty() );
}

BregmanProjectionFunctional::BregmanProjectionFunctional(
		const InverseProblem_ptr_t &_problem,
		const std::vector<SpaceElement_ptr_t> &_U,
		const std::vector<double> &_alpha
		) :
	dualpower(_problem->SourceSpace->getDualSpace()->getDualityMapping()->getPower()),
	dualnorm(*_problem->SourceSpace->getDualSpace()->getNorm()),
	J_q(*_problem->SourceSpace->getDualSpace()->getDualityMapping()),
	U(_U),
	alpha(_alpha),
	resx((*U.begin())->getSpace()->createElement()),
	dual_resx(resx->getSpace()->getDualSpace()->createElement()),
	zeroVec((*U.begin())->getSpace()->createElement())
{
	assert ( U.size() == alpha.size() );
	assert( !U.empty() );
}

void BregmanProjectionFunctional::updateDualIterate(
		const std::vector<double> &_t,
		const SpaceElement_ptr_t &_dualx) const
{
	resx->setZero();
	for (size_t dim = 0; dim < _t.size(); ++dim)
		resx->scaledAddition(-1.*_t[dim], U[dim]);
	*resx += _dualx;
}

double BregmanProjectionFunctional::operator()(
		const std::vector<double> &_t,
		const SpaceElement_ptr_t &_dualx
		) const
{
	assert ( _t.size() == U.size() );
	assert( (*U.begin())->getSpace() == _dualx->getSpace() );
	// x=x-U*t;
	updateDualIterate(_t, _dualx);
	// fval=1/q*norm(x,p)^q+alpha'*t;
	double alpha_times_t = 0.;
	alpha_times_t = std::inner_product(
			_t.begin(), _t.end(),
			alpha.begin(),
			alpha_times_t);
	const double fval =
			1./dualpower * ::pow(dualnorm(resx), dualpower)
			+ alpha_times_t;
	if (isnan(fval)) {
		LOG(error, "Encountered NaN in function calculation with " << *_dualx
//				<< " at " << _t << " and argument ";
				);
	}
	return fval;
}

std::vector<double> BregmanProjectionFunctional::gradient(
		const std::vector<double> &_t,
		const SpaceElement_ptr_t &_dualx
		) const
{
	assert ( _t.size() == U.size() );
	assert( (*U.begin())->getSpace() == _dualx->getSpace() );
	// x=x-U*t;
	updateDualIterate(_t, _dualx);
	J_q(resx, dual_resx);
	std::vector<double> gval(alpha);
	for (size_t i=0;i<_t.size();++i) {
		gval[i] -= U[i] * dual_resx;
		if (isnan(gval[i])) {
			LOG(error, "Encountered NaN in gradient calculation  with " << *_dualx
//					<< " at " << _t << " and argument ";
					);
		}
	}
	return gval;
}
