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
 * NemirovskyDirection.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NemirovskyDirection.hpp"

#include <cassert>

#include "Log/Logging.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"

NemirovskyDirection::NemirovskyDirection(
		const NormedSpace_ptr_t &_SearchDirectionSpace_ptr) :
		Searchspace(_SearchDirectionSpace_ptr, 2)
{}

void NemirovskyDirection::update(
		const SpaceElement_ptr_t &_newdir,
		const double _alpha,
		const SpaceElement_ptr_t &_dual_iterate,
		const SpaceElement_ptr_t &_iterate
		)
{
	assert( _newdir->getSpace() == SearchDirectionSpace_ptr );
	assert( _dual_iterate->getSpace() == SearchDirectionSpace_ptr );

	// update search direction
	*(U[0]) = _newdir;
	alphas[0] = _alpha;

	// update Nemirovsky direction
	*(U[1]) = _dual_iterate;
	alphas[1] = _dual_iterate * _iterate;
}
