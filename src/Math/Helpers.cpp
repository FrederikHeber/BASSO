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

double Helpers::sign(const double _x)
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


