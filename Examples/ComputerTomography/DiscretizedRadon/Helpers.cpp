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
 * Helpers.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "ComputerTomography/DiscretizedRadon/Helpers.hpp"

#include <algorithm>

namespace detail {

	std::vector<double> calculatePixelCenter(
			const double _h,
			const unsigned int _num_pixel)
	{
		std::vector<double> values(_num_pixel);
		for (unsigned int pixel_x = 0; pixel_x < _num_pixel; ++pixel_x)
			values[pixel_x] = -1.+(pixel_x+.5)*_h;
		return values;
	}

	std::vector< point_t > calculateOmegas(const unsigned int _num_angles)
	{
		std::vector< point_t > omegas(_num_angles);
		const double deltaphi = M_PI/_num_angles;
		std::vector<double> angles(_num_angles);
		std::generate(
				angles.begin(), angles.end(),
				AngleIncrementer(deltaphi));
		std::transform(
				angles.begin(), angles.end(),
				omegas.begin(),
				calculateOmega());
		return omegas;
	}

	point_t calculateOmega::operator()(const double _phi) const
	{
		point_t omega;
		omega[0] = cos(_phi);
		omega[1] = sin(_phi);
		return omega;
	}

	AngleIncrementer::AngleIncrementer(const double _delta) :
		angle(-_delta),
		delta(_delta)
	{}
	double AngleIncrementer::operator()()
	{
		angle += delta;
		return angle;
	}

}; /* namespace detail */
