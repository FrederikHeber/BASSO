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
