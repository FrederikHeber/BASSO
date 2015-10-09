/*
 * MatrixToPNG.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <boost/chrono.hpp>
#include <png++/png.hpp>

#include "Options/MatrixToPNGOptions.hpp"

#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"

int main (int argc, char *argv[])
{
	// starting timing
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	MatrixToPNGOptions opts;
	opts.init();

	// parse options
	opts.parse(argc, argv);

	if (opts.showHelpConditions(argv[0]))
		return 1;

	// set verbosity level
	opts.setVerbosity();

	if (!opts.checkSensibility())
		return 255;
	opts.setSecondaryValues();

	// parse vector files into instances
	Eigen::VectorXd matrix;
	const double num_pixels = opts.num_pixel_x * opts.num_pixel_y;
	{
		using namespace MatrixIO;

		{
			std::ifstream ist(opts.matrix_file.string().c_str());
			if (ist.good())
				try {
					ist >> matrix;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse rhs from " << opts.matrix_file.string() << std::endl;
					return 255;
				}
			else {
				std::cerr << "Failed to open " << opts.matrix_file.string() << std::endl;
				return 255;
			}
		}
	}
	assert( matrix.rows() == num_pixels);

	if (!boost::filesystem::exists(opts.image_file)) {
		// find bounds
		const double min = matrix.minCoeff();
		const double max = matrix.maxCoeff();
		const double length = fabs(max - min);
		png::image< png::rgb_pixel > image(opts.num_pixel_x, opts.num_pixel_y);
		for (size_t y = 0; y < image.get_height(); ++y)
		{
		 for (size_t x = 0; x < image.get_width(); ++x)
		 {
			 unsigned int i;
			 unsigned int j;
			 unsigned int multiplier = image.get_width();

			 if (opts.LeftToRight)
				 i = x;
			 else
				 i = image.get_width() - x - 1;

			 if (opts.BottomToTop)
				 j = y;
			 else
				 j = image.get_height() - y - 1;

			 if (opts.Flip) {
				 multiplier = image.get_height();
				 std::swap(i,j);
			 }

			 const unsigned int value = 255.*(matrix[i+j*multiplier] - min)/length;
			 image[y][x] = png::rgb_pixel(value, value, value);
			 // non-checking equivalent of image.set_pixel(x, y, ...);
		 }
		}
		image.write(opts.image_file.string());
	} else {
		BOOST_LOG_TRIVIAL(error) << "The given output file exists.";
	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	BOOST_LOG_TRIVIAL(info) << "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start)
			<< ".";

	// exit
	return 0;
}



