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
 * SingularValueDecomposition_impl.cpp
 *
 *  Created on: Oct 5, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SingularValueDecomposition_impl.hpp"

SingularValueDecomposition_impl::SingularValueDecomposition_impl(
		const Eigen::MatrixXd &_matrix) :
	svd(_matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV))
{}

