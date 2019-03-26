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
 * LinearDependencyChecker.cpp
 *
 *  Created on: Feb 10, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "LinearDependencyChecker.hpp"

#include <Eigen/Dense>

#include "Log/Logging.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

bool LinearDependencyChecker::operator()(
		const vectors_t &_vectors)
{
	if (_vectors.empty()) {
		LOG(warning, "No vectors given whose linear dependency to check.");
		return true;
	}
	if (_vectors.size() == 1) {
		const bool result = _vectors[0]->isZero(BASSOTOLERANCE);
		LOG(debug, "Checking Linear Dependency: Rank is "
				<< (result ? "0" : "1") << " with " << _vectors.size() << " vectors.");
		return result;
	}
	Eigen::MatrixXd matrix(
			_vectors[0]->getSpace()->getDimension(),
			_vectors.size());
	for (size_t index = 0; index < _vectors.size(); ++index)
		matrix.col(index) =
				RepresentationAdvocate::get(_vectors[index]);
	const unsigned int rank = matrix.colPivHouseholderQr().rank();
	LOG(debug, "Checking Linear Dependency: Rank is " << rank << " with " << _vectors.size() << " vectors.");
	return rank != _vectors.size();
}


