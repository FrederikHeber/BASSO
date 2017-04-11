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
 * InverseProblem.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InverseProblem.hpp"

#include <cassert>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

InverseProblem::InverseProblem(
		const Mapping_ptr_t &_A,
		const NormedSpace_ptr_t &_SourceSpace,
		const NormedSpace_ptr_t &_TargetSpace,
		const SpaceElement_ptr_t &_y
		) :
	SourceSpace(_SourceSpace),
	DualSourceSpace(_SourceSpace->getDualSpace()),
	TargetSpace(_TargetSpace),
	DualTargetSpace(_TargetSpace->getDualSpace()),
	A(_A),
	A_t(_A->getAdjointMapping()),
	x(_SourceSpace->createElement()),
	y(_y)
{
	assert( A->getTargetSpace() == y->getSpace() );
}

