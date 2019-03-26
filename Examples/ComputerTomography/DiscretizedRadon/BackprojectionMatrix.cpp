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
 * BackprojectionMatrix.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: heber
 */


#include <ComputerTomography/DiscretizedRadon/BackprojectionMatrix.hpp>
#include "BassoConfig.h"

#include <fstream>
#include <vector>

#include "ComputerTomography/DiscretizedRadon/Helpers.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"


BackprojectionMatrix::BackprojectionMatrix(
		const unsigned int _num_pixel_x,
		const unsigned int _num_pixel_y,
		const unsigned int _num_angles,
		const unsigned int _num_offsets) :
		num_pixel_x(_num_pixel_x),
		num_pixel_y(_num_pixel_y),
		num_angles(_num_angles),
		num_offsets(_num_offsets),
		matrix(_num_pixel_x*_num_pixel_y, _num_angles*_num_offsets)
{
	matrix.setZero();

	// precalculate omegas and pixel centers
	std::vector< point_t > omegas = detail::calculateOmegas(num_angles);
	const double h_x = 2./(double)num_pixel_x;
	const double h_y = 2./(double)num_pixel_y;
	std::vector< double > pixelcenter_x =
			detail::calculatePixelCenter(h_x, num_pixel_x);
	std::vector< double > pixelcenter_y =
			detail::calculatePixelCenter(h_y, num_pixel_y);

	// go through every pixel
	assert( _num_offsets % 2 == 1);
	const int half_offsets =
			_num_offsets > 1 ?
					(_num_offsets-1) / 2 : 0;
	const double delta = 1./(double)half_offsets;
	const double prefactor = M_PI/(double)num_angles;
	for (unsigned int pixel_x = 0; pixel_x < num_pixel_x; ++pixel_x) {
		for (unsigned int pixel_y = 0; pixel_y < num_pixel_y; ++pixel_y) {
			const unsigned int row_index = pixel_y + (pixel_x*num_pixel_y);
			const point_t x(pixelcenter_x[pixel_x], pixelcenter_y[pixel_y]);
//			if (x.squaredNorm() - 1. > BASSOTOLERANCE)
//				continue;
			for (unsigned int angle = 0; angle < _num_angles; ++angle) {
				const double s =
						std::min(1.,std::max(-1.,omegas[angle].dot(x)));
				const double t = s/delta;
				const int k = floor(t);
				const double u = t - (double)k;
				assert( (u > -BASSOTOLERANCE) && (u-1. < BASSOTOLERANCE) );
				const unsigned int col_index =
						(k+half_offsets) + (angle*num_offsets);
				const unsigned int next_col_index =
						col_index + 1;
				if ((k+half_offsets) >= 0)
					matrix.coeffRef( row_index, col_index) += prefactor*(1.-u);
				if ((k+half_offsets+1) < (int)num_offsets)
					matrix.coeffRef( row_index, next_col_index) += prefactor*u;
			}
		}
	}
}

int BackprojectionMatrix::load(const boost::filesystem::path &_file)
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
