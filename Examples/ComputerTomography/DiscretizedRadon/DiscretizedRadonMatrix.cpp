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
 * DiscretizedRadonMatrix.cpp
 *
 *  Created on: Jul 6, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "DiscretizedRadonMatrix.hpp"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iterator>
#include <iostream>
#include <numeric>
#include <sstream>

#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"

DiscretizedRadonMatrix::DiscretizedRadonMatrix(
		const unsigned int _num_pixel_x,
		const unsigned int _num_pixel_y,
		const unsigned int _num_angles,
		const unsigned int _num_offsets) :
		num_pixel_x(_num_pixel_x),
		num_pixel_y(_num_pixel_y),
		num_angles(_num_angles),
		num_offsets(_num_offsets),
		matrix(_num_angles*_num_offsets, _num_pixel_x*_num_pixel_y),
		debugflag(false)
{
	matrix.setZero();

//	const double h = 2./(double)_num_offsets;
	// all parameters must be positive
	assert( _num_pixel_x > 0);
	assert( _num_pixel_y > 0);
	assert( _num_angles > 0);
	assert( _num_offsets > 0);
	// num_offsets must be odd
	assert( _num_offsets % 2 == 1);
	const int half_offsets =
			_num_offsets > 1 ?
					(_num_offsets-1) / 2 : 0;
	if (debugflag) {
		LOG(debug, "half_offsets = " << half_offsets);
	}
	const double q =
			(half_offsets != 0) ?
					1./(double)half_offsets : 2.;
	if (debugflag) {
		LOG(debug, "deltas = " << q);
	}
	const double deltaphi = M_PI/_num_angles;
	if (debugflag) {
		LOG(debug, "deltaphi = " << deltaphi*180./M_PI);
	}

	// a_jm is the length of the intersection of line j with pixel m
	for (unsigned int angle = 0; angle < _num_angles; ++angle) {
		const double phi = angle*deltaphi;
		if (debugflag) {
			LOG(debug, "phi(" << angle << ") = " << phi);
		}
		for (int offset = -half_offsets; offset <= half_offsets; ++offset) {
			const double s = q * offset;
			if (debugflag) {
				LOG(debug, "s(" << offset+half_offsets << ") = " << s);
			}

			// get all intersections points for this line
			intersections_t intersections =
					calculatePixelBoundaryIntersections(phi,s);
			if (debugflag)
			{
				std::stringstream output;
				std::copy(intersections.begin(), intersections.end(),
						std::ostream_iterator<point_t>(output, "\n"));
				LOG(debug, "intersections = \n" << output.str());
			}

			// remove identical and illegal points
			removeIdenticalAdjacentPoints(intersections);
//			removeIllegalPoints(intersections);

			if (intersections.size() > 1) {
				// convert to pixels
				pixels_t pixels = PointsToPixels(intersections);
				if (debugflag)
				{
					std::stringstream output;
					for (pixels_t::const_iterator pixeliter = pixels.begin();
							pixeliter != pixels.end(); ++pixeliter) {
						const unsigned int col_index =
								(*pixeliter)[1] + ((*pixeliter)[0] * num_pixel_y);
						output << col_index+1 << "\n";
					}
					LOG(debug, "pixels = \n" << output.str());
				}

				// flip pixels on upper/right image edges
				correctBoundaryPixels(angle, offset, pixels);
				if (debugflag)
				{
					std::stringstream output;
					for (pixels_t::const_iterator pixeliter = pixels.begin();
							pixeliter != pixels.end(); ++pixeliter) {
						const unsigned int col_index =
								(*pixeliter)[1] + ((*pixeliter)[0] * num_pixel_y);
						output << col_index+1 << "\n";
					}
					LOG(debug, "pixels = \n" << output.str());
				}

				// calculate lengths and add onto pixels
				calculateLengths(angle, offset + half_offsets, pixels, intersections);
			}
		}
	}
}

int DiscretizedRadonMatrix::load(const boost::filesystem::path &_file)
{
	using namespace MatrixIO;

	{
		std::ifstream ist(_file.string().c_str());
		if (ist.good())
			try {
				ist >> matrix;
			} catch (MatrixIOStreamEnded_exception &e) {
				std::cerr << "Failed to fully parse matrix from " << _file.string() << std::endl;
				return 255;
			}
		else {
			std::cerr << "Failed to open " << _file.string() << std::endl;
			return 255;
		}
	}

	return 0;
}

/** Line is parametrized as
 * \f$ (x \cdot -\sin(\phi) + s \cdot \cos(\phi) , x \cdot \cos(\phi) + s \sin(\phi)) \f$
 *
 * @param _phi angle of the line with respect to x axis
 * @param _s offset of the line in orthogonal direction to angle vector
 * @param _x length parameter of the line
 * @return coordinates in dim-dimensional space
 */
static
point_t convertLineToPoint(
		const double _cphi,
		const double _sphi,
		const double _s,
		const double _x)
{
	point_t point;
	assert( point.size() == 2 );
	point[0] = _x * (-_sphi) + _s * _cphi;
	point[1] = _x * _cphi + _s * _sphi;

	return point;
}

void DiscretizedRadonMatrix::calculateLengths(
		unsigned int _angle,
		unsigned int _offset,
		const pixels_t &_pixels,
		const intersections_t &_intersections)
{
	const unsigned int row_index =
			_offset + (_angle * num_offsets);
	const double totallength =
			(*_intersections.begin() - *_intersections.rbegin()).norm();
	pixels_t::const_iterator pixeliter = _pixels.begin();
	intersections_t::const_iterator advanceiter = _intersections.begin();
	intersections_t::const_iterator pointiter = advanceiter++;
	double length_sum = 0.;
	for (; advanceiter != _intersections.end();
			pointiter = advanceiter++, ++pixeliter) {
		const double length = ((*pointiter) - (*advanceiter)).norm();
		const unsigned int col_index =
				(*pixeliter)[1] + ((*pixeliter)[0] * num_pixel_y);
		matrix.coeffRef( row_index, col_index ) = length;
		length_sum += length;
		if (debugflag) {
			LOG(debug, "a_{" << row_index << "," << col_index << "} = " << length);
		}
	}
	assert( fabs(length_sum - totallength) < BASSOTOLERANCE );
	assert( pixeliter == _pixels.end() );
}

DiscretizedRadonMatrix::intersections_t
DiscretizedRadonMatrix::calculatePixelBoundaryIntersections(
		const double _phi,
		const double _s) const
{
	intersections_t intersections;

	const double cphi = cos(_phi);
	const double sphi = sin(_phi);
	const double h_x = 2./(double)num_pixel_x;
	const double h_y = 2./(double)num_pixel_y;
	if (fabs(sphi) > BASSOTOLERANCE)
		// go one further to also calculate the end point of the last pixel
		for (unsigned int pixel_x = 0; pixel_x <= num_pixel_x; ++pixel_x) {
			/// calculate x coordinate of intersection with vertical line
			const double pos_x = -1. + h_x*pixel_x;

			/// calculate intersection
			const double length_pos_x = (-pos_x + _s * cphi) / sphi;

			/// add to sorted array. This makes sense because points
			/// all are on a line, i.e. 1d manifold, and hence can be
			/// aligned easily and sensibly.
			const point_t point_x =
					convertLineToPoint(cphi, sphi, _s, length_pos_x);
			/// only add if within [-1,-1]^dim domain
			if (point_x.lpNorm<Eigen::Infinity>() <= 1.)
				intersections.insert(point_x);
		}
	if (fabs(cphi) > BASSOTOLERANCE)
		// go one further to also calculate the end point of the last pixel
		for (unsigned int pixel_y = 0; pixel_y <= num_pixel_y; ++pixel_y) {
			/// calculate y coordinate of intersection with horizontal line
			const double pos_y = -1. + h_y*pixel_y;

			/// calculate intersections
			const double length_pos_y = (pos_y - _s * sphi) / cphi;

			/// add to sorted array. This makes sense because points
			/// all are on a line, i.e. 1d manifold, and hence can be
			/// aligned easily and sensibly.
			const point_t point_y =
					convertLineToPoint(cphi, sphi, _s, length_pos_y);
			/// only add if within [-1,-1]^dim domain
			if (point_y.lpNorm<Eigen::Infinity>() <= 1.)
				intersections.insert(point_y);
		}
	if ((0) && (intersections.size() > 0)) {
		LOG(debug, "Line starts at " << _s << " and slants with (" << cphi << "," << sphi << ").");
		std::stringstream output;
		std::copy(intersections.begin(), intersections.end(),
				std::ostream_iterator<point_t>(output, "\n"));
		LOG(debug, "We have the following intersections at Pixel from : " << std::endl << output.str());
	}
	return intersections;
}

void DiscretizedRadonMatrix::removeIdenticalAdjacentPoints(
		intersections_t &_intersections) const
{
	if (!_intersections.empty()) {
		intersections_t::iterator advanceiter = _intersections.begin();
		intersections_t::iterator pointiter = advanceiter++;
		for (; advanceiter != _intersections.end();) {
			const double length = ((*pointiter) - (*advanceiter)).squaredNorm();
			if (length < BASSOTOLERANCE) {
				_intersections.erase(advanceiter);
				advanceiter = pointiter;
				if (++advanceiter == _intersections.end())
					break;
			} else {
				pointiter = advanceiter++;
			}
		}
	}
}

void DiscretizedRadonMatrix::removeIllegalPoints(
		intersections_t &_intersections) const
{
	if (!_intersections.empty()) {
		intersections_t::iterator advanceiter = _intersections.begin();
		intersections_t::iterator pointiter = advanceiter++;
		for (; advanceiter != _intersections.end();
				pointiter = advanceiter++) {
			if (((*pointiter)[0]+BASSOTOLERANCE >= 1.) || ((*pointiter)[1]+BASSOTOLERANCE >= 1.))
				_intersections.erase(pointiter);
		}
		// also check last element
		if (((*pointiter)[0]+BASSOTOLERANCE >= 1.) || ((*pointiter)[1]+BASSOTOLERANCE >= 1.))
			_intersections.erase(pointiter);
	}
}

void DiscretizedRadonMatrix::correctBoundaryPixels(
		const unsigned int _angle,
		const unsigned int _offset,
		pixels_t &_pixels) const
{
	const double deltaphi = M_PI/num_angles;
	const double phi = _angle*deltaphi;
	if ((fabs(sin(phi)) < BASSOTOLERANCE) && (_offset == (num_offsets-1)/2))
		for(pixels_t::iterator pixeliter = _pixels.begin();
			pixeliter != _pixels.end(); ++pixeliter)
			if ((*pixeliter)[0] == num_pixel_x)
				(*pixeliter)[0] -= 1;
	if ((fabs(cos(phi)) < BASSOTOLERANCE) && (_offset == (num_offsets-1)/2))
		for(pixels_t::iterator pixeliter = _pixels.begin();
			pixeliter != _pixels.end(); ++pixeliter)
			if ((*pixeliter)[1] == num_pixel_y)
				(*pixeliter)[1] -= 1;
}

DiscretizedRadonMatrix::pixels_t DiscretizedRadonMatrix::PointsToPixels(
		const intersections_t &_intersections) const
{
	pixel_t zeropixel;
	zeropixel[0]=-1;
	zeropixel[1]=-1;
	const double h_x = 2./(double)num_pixel_x;
	const double h_y = 2./(double)num_pixel_y;
	pixels_t pixels(_intersections.size()-1, zeropixel);
	pixels_t::iterator pixeliter = pixels.begin();
	intersections_t::const_iterator advanceiter = _intersections.begin();
	intersections_t::const_iterator pointiter = advanceiter++;
	for (; advanceiter != _intersections.end();
			pointiter = advanceiter++, ++pixeliter) {
		// center between adjacent intersections is always
		point_t center = (*pointiter);
		center += (*advanceiter);
		center *= .5;
		(*pixeliter)[0] = (unsigned int)floor((1.+center[0])/h_x);
		(*pixeliter)[1] = (unsigned int)floor((1.+center[1])/h_y);
//		if ((*pixeliter)[0] >= num_pixel_x)
//			(*pixeliter)[0] = num_pixel_x-1;
//		if ((*pixeliter)[1] >= num_pixel_y)
//			(*pixeliter)[1] = num_pixel_y-1;
	}

	return pixels;
}

std::ostream& operator<<(
		std::ostream &ost, const DiscretizedRadonMatrix::pixel_t &_pixel)
{
	ost << "[" << _pixel[0] << "," << _pixel[1] << "]";
	return ost;
}
