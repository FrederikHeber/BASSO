/*
 * Helpers.cpp
 *
 *  Created on: Jul 22, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "ComputerTomographyBase/DiscretizedRadon/Helpers.hpp"

#include <algorithm>

namespace detail {

	point_t getPixelCenter(
			const unsigned int _i,
			const unsigned int _j,
			const double _h_x,
			const double _h_y
			)
	{
		point_t point;
		point[0] = -1.+(_i+.5)*_h_x;
		point[1] = -1.+(_j+.5)*_h_y;
		return point;
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
