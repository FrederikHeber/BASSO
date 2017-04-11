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
 * SingularValueDecomposition.cpp
 *
 *  Created on: Oct 5, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <cassert>

#include <Eigen/Dense>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/SingularValueDecomposition.hpp"
#include "Minimizations/Mappings/SingularValueDecomposition_impl.hpp"
#include "Minimizations/Norms/Norm.hpp"

SpaceElement_ptr_t
SingularValueDecomposition::solve(SpaceElement_ptr_t _rhs) const
{
	// SVD only gives true solution for l2 norm
	assert( _rhs->getSpace()->getNorm()->getPvalue() == 2. );
	const Eigen::VectorXd &rhsvector =
			RepresentationAdvocate::get(_rhs);
	const Eigen::VectorXd solution = pimpl->svd.solve(rhsvector);
	return ElementCreator::create(_rhs->getSpace(), solution);
}
