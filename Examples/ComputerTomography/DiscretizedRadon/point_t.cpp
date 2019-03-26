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
 * point_t.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "point_t.hpp"

#include <iostream>

bool point_t::operator<(
		const point_t &_other) const
{
	for (unsigned int i=0;i<dim;++i) {
		if ((*this)[i] > _other[i])
			return false;
		else if ((*this)[i] < _other[i])
			return true;
	}
	return false;
}

point_t::point_t(const double _x, const double _y) :
	 Eigen::Vector2d(_x,_y)
{}

point_t::point_t(
		const Eigen::Vector2d &_p) :
	 Eigen::Vector2d(_p)
{}

point_t::point_t(
		const Eigen::Transpose<Eigen::Vector2d> &_p) :
	 Eigen::Vector2d(_p)
{}

point_t::point_t()
{}

std::ostream& operator<<(
		std::ostream &_ost, const point_t &_point)
{
	_ost << "[" << _point[0] << "," << _point[1] << "]";
	return _ost;
}
