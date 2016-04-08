/*
 * MatrixToPNG.cpp
 *
 *  Created on: Oct 9, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <fstream>

#include <boost/chrono.hpp>
#include <png++/png.hpp>

#include "Options/MatrixToPNGOptions.hpp"

#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"

int main (int argc, char *argv[])
{
	// show program information
	showVersion(std::string(argv[0]));
	showCopyright();

	MatrixToPNGOptions opts;
	opts.init();

	// parse options
	try {
		opts.parse(argc, argv);
	} catch (std::exception &e) {
		std::cerr << "An error occurred: "
				<< e.what()
				<< std::endl;
		return 255;
	}


	if (opts.showHelpConditions(argv[0]))
		return 1;

	// set verbosity level
	opts.setVerbosity();

	if (!opts.checkSensibility())
		return 255;
	opts.setSecondaryValues();

	// starting timing
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	// parse vector files into instances
	Eigen::VectorXd matrix;
	const double num_pixels = opts.num_pixel_x * opts.num_pixel_y;
	if (!MatrixIO::parse(opts.matrix_file.string(), "matrix", matrix))
		return 255;
	assert( matrix.rows() == num_pixels);

	if (!boost::filesystem::exists(opts.image_file)) {
		// find bounds
		const double min = matrix.minCoeff();
		const double max = matrix.maxCoeff();
		const double length = fabs(max - min);
		BOOST_LOG_TRIVIAL(info)
				<< "Data range is [" << min << ":" << max << "]";

		png::image< png::rgb_pixel > *image = NULL;
		switch(opts.Rotate) {
		case 0:
		case 2:
			image = new png::image< png::rgb_pixel >(opts.num_pixel_y, opts.num_pixel_x);
			break;
		case 1:
		case 3:
			image = new png::image< png::rgb_pixel >(opts.num_pixel_x, opts.num_pixel_y);
			break;
		default:
			assert(0); /* Case cannot happen */
			break;
		}
		assert( image != NULL );
		unsigned int multiplier = !opts.Flip ? opts.num_pixel_x : opts.num_pixel_y;
		for (size_t y = 0; y < opts.num_pixel_y; ++y)
		{
		 for (size_t x = 0; x < opts.num_pixel_x; ++x)
		 {
			 unsigned int i;
			 unsigned int j;

			 if (opts.LeftToRight)
				 i = x;
			 else
				 i = opts.num_pixel_x - x - 1;

			 if (opts.BottomToTop)
				 j = y;
			 else
				 j = opts.num_pixel_y - y - 1;

			 if (opts.Flip)
				 std::swap(i,j);

			 const unsigned int value = 255.*(matrix[i+j*multiplier] - min)/length;
			 const png::rgb_pixel pixel_value = png::rgb_pixel(value, value, value);
			 switch(opts.Rotate) {
			 case 0:
				 (*image)[x][y] = pixel_value;
				 break;
			 case 1:
				 (*image)[opts.num_pixel_y-1-y][x] = pixel_value;
				 break;
			 case 2:
				 (*image)[opts.num_pixel_x-1-x][opts.num_pixel_y-1-y] = pixel_value;
				 break;
			 case 3:
				 (*image)[y][opts.num_pixel_x-1-x] = pixel_value;
				 break;
			 default:
				 assert(0); /* case should not happen */
				 break;
			 }
			 // non-checking equivalent of image->set_pixel(x, y, ...);
		 }
		}
		std::ofstream output(opts.image_file.string().c_str());
		image->write_stream(output);
		output.close();
		if (image != NULL)
			delete image;
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



