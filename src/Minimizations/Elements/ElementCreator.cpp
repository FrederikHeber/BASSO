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
 * ElementCreator.cpp
 *
 *  Created on: Dec 11, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "ElementCreator.hpp"

const SpaceElement_ptr_t ElementCreator::create(
		const NormedSpace_weakptr_t _space,
		const Eigen::VectorXd &_vector)
{
	return create(*NormedSpace_ptr_t(_space), _vector);
}

const SpaceElement_ptr_t ElementCreator::create(
		const NormedSpace &_space,
		const Eigen::VectorXd &_vector)
{
	SpaceElement_ptr_t returnelement = _space.createElement();
	*returnelement = _vector;
	return returnelement;
}
