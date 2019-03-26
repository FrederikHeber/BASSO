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
 * DynamicRegularizedL1NormStepWidth.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "DynamicRegularizedL1NormStepWidth.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/NormFactory.hpp"

DynamicRegularizedL1NormStepWidth::DynamicRegularizedL1NormStepWidth(
		const InverseProblem_ptr_t &_problem,
		const Mapping_ptr_t &_J_r) :
	problem(_problem),
	J_r(_J_r),
	A_adjoint(_problem->A->getAdjointMapping()),
	l2norm_Y(NormFactory::create(
			"lp",
			problem->A->getTargetSpace(),
			NormFactory::args_t(1, boost::any(2.)))),
	l2norm_DualX(NormFactory::create(
			"dual_lp",
			problem->x->getSpace()->getDualSpace(),
			NormFactory::args_t(1, boost::any(2.))))
{}

const double DynamicRegularizedL1NormStepWidth::operator()(
		const SpaceElement_ptr_t &_dualx,
		const SpaceElement_ptr_t &_u,
		const SpaceElement_ptr_t &_solution,
		const SpaceElement_ptr_t &_residual,
		const double _residuum,
		const double _TolX,
		const double _alpha
		) const
{
	const LinearMapping &A = static_cast<const LinearMapping &>(*problem->A);
	const LinearMapping &A_t = static_cast<const LinearMapping &>(*A_adjoint);
	const SpaceElement_ptr_t upper = A * _solution - problem->y;
	const SpaceElement_ptr_t lower = A_t * (*J_r)(upper);
	const double numerator = (*l2norm_Y)(upper);
	const double denominator = (*l2norm_DualX)(lower);
	if (fabs(denominator) > BASSOTOLERANCE) {
		const double norm_value = numerator/denominator;
		return norm_value*norm_value;
	} else {
		return 0.;
	}
}

