/*
 * BackprojectionMatrix.cpp
 *
 *  Created on: Jul 8, 2015
 *      Author: heber
 */


#include <ComputerTomographyBase/DiscretizedRadon/BackprojectionMatrix.hpp>
#include "BassoConfig.h"

#include <algorithm>
#include <fstream>
#include <vector>

#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"

/** Returns the center of the pixel indicated by the indices \a _i
 * and \a _j, this is times 1 over h_x/h_y.
 *
 * @param _i index for x
 * @param _j index for y
 * @param _h_x length of a pixel in x direction
 * @param _h_y length of a pixel in y direction
 * @return pixel center coordinates
 */
static point_t getPixelCenter(
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

struct calculateOmega
{
	point_t operator()(const double _phi) const
	{
		point_t omega;
		omega[0] = cos(_phi);
		omega[1] = sin(_phi);
		return omega;
	}
};

struct AngleIncrementer
{
	AngleIncrementer(const double _delta) :
		angle(-_delta),
		delta(_delta)
	{}
	double operator()()
	{
		angle += delta;
		return angle;
	}

private:
	double angle;
	double delta;
};

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

	// precalculate omegas
	std::vector<double> angles(num_angles);
	std::generate(
			angles.begin(), angles.end(),
			AngleIncrementer(M_PI/(double)num_angles));
	std::vector< point_t > omegas(num_angles);
	std::transform(
			angles.begin(), angles.end(),
			omegas.begin(),
			calculateOmega());

	// go through every pixel
	assert( _num_offsets % 2 == 1);
	const int half_offsets = (_num_offsets-1) / 2;
	const double delta = 1./(double)half_offsets;
	const double h_x = 2./(double)num_pixel_x;
	const double h_y = 2./(double)num_pixel_y;
	const double prefactor = M_PI/(double)num_angles;
	for (unsigned int pixel_x = 0; pixel_x < num_pixel_x; ++pixel_x) {
		for (unsigned int pixel_y = 0; pixel_y < num_pixel_y; ++pixel_y) {
			const unsigned int row_index = pixel_y + (pixel_x*num_pixel_y);
			const point_t x = getPixelCenter(pixel_x, pixel_y, h_x, h_y);
//			if (x.squaredNorm() - 1. > BASSOTOLERANCE)
//				continue;
			for (unsigned int angle = 0; angle < _num_angles; ++angle) {
				const double s = omegas[angle].dot(x);
				const double t = s/delta;
				const int k = floor(t);
				const double u = t - (double)k;
				assert( (u > -BASSOTOLERANCE) && (u-1. < BASSOTOLERANCE) );
				const unsigned int col_index =
						(k+half_offsets) + (angle*num_offsets);
				if (((k+half_offsets) < 0)
						|| ((k+half_offsets+1) >= (int)num_offsets))
					continue;
				const unsigned int next_col_index =
						col_index + 1;
				matrix( row_index, col_index) += prefactor*(1.-u);
				matrix( row_index, next_col_index) += prefactor*u;
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
