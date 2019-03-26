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
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "Helpers.hpp"
#include "MathExceptions.hpp"

#include <cmath>

Eigen::VectorXd Helpers::circshift(const Eigen::VectorXd &_x, const int shift)
{
	if (shift == 0)
	return _x;
	else {
		// for the moment we just copy the entries and do not use any ...
		// TODO: fancy mem copying of blocks
		const int size = _x.innerSize();
		if (abs(shift) > size)
			throw MathIllegalValue_Error()
				<< MathIllegalValue_name("shift");
		Eigen::VectorXd shiftedx(size);
		for (int i=0;i<size; ++i) {
			// add size to prevent negative numbers
			const int index = (i+shift+size)%size;
			shiftedx[index] = _x[i];
		}
		return shiftedx;
	}
}

double Helpers::sign(const double& _x)
{
	if (fabs(_x) < BASSOTOLERANCE)
		return 0.;
	else
		return (_x/fabs(_x));
}

Eigen::VectorXd Helpers::signum(const Eigen::VectorXd &_x)
{
	// signbit returns 1 (negative) or 0 (positive)
	const int size = _x.innerSize();
	Eigen::VectorXd signedx(size);
	for (int i=0;i<size; ++i) {
		signedx[i] = Helpers::sign(_x[i]);
	}
	return signedx;
}


