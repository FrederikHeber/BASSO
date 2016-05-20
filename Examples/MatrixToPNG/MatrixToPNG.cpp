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

#include "Colortables/ColorTable.hpp"
#include "Options/MatrixToPNGOptions.hpp"

#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"

png::rgb_pixel getPixelValue(
		const double _value,
		const double _min,
		const double _max,
		const Eigen::MatrixXd &_colortable)
{
	png::rgb_pixel pixel_value;
	if ((_colortable.rows() % 2 != 0) && (_min < 0)) {
		// odd number of rows: middle most row gives "zero" color
		const int no_sections = floor(_colortable.rows()/2);
		const double length = (_value > 0) ? fabs(_max) : fabs(_min);
		double temp = 0.;
		if (_value > 0)
			temp = (double)no_sections*_value/length;
		else
			temp = (double)no_sections*fabs(_value - _min)/length;
		temp = std::max(std::min(temp, (double)no_sections), 0.);
		const int section = floor(temp);
		const double fraction = temp - section;
		int offset = 0;
		if (_value > 0) { // use second half of color table
			offset = no_sections;
		}
		const int either_section = offset + section;
		int or_section = offset + section + 1;
		if (or_section == _colortable.rows())
			or_section = either_section;
		assert( (either_section >=0 ) && (either_section < _colortable.rows()));
		assert( (or_section >=0 ) && (or_section < _colortable.rows()));
		pixel_value = png::rgb_pixel(
				255*(_colortable(either_section,0)*(1.-fraction) + _colortable(or_section,0)*fraction),
				255*(_colortable(either_section,1)*(1.-fraction) + _colortable(or_section,1)*fraction),
				255*(_colortable(either_section,2)*(1.-fraction) + _colortable(or_section,2)*fraction)
				);
	} else {
		const int no_sections = _colortable.rows()-1;
		const double length = fabs(_max - _min);
		const double temp = (double)no_sections*fabs(_value - _min)/length;
		const int section = floor(temp);
		const double fraction = temp - section;
		const int either_section = section;
		int or_section = section + 1;
		if (or_section == _colortable.rows())
			or_section = either_section;
		assert( (either_section >=0 ) && (either_section < _colortable.rows()));
		assert( (or_section >=0 ) && (or_section < _colortable.rows()));
		pixel_value = png::rgb_pixel(
				255*(_colortable(either_section,0)*(1.-fraction) + _colortable(or_section,0)*fraction),
				255*(_colortable(either_section,1)*(1.-fraction) + _colortable(or_section,1)*fraction),
				255*(_colortable(either_section,2)*(1.-fraction) + _colortable(or_section,2)*fraction)
				);
	}
	return pixel_value;
}

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

	// get color table
	ColorTable table;
	Eigen::MatrixXd colortable;
	if (!opts.Colorize.empty()) {
		if (table.isKeyPresent(opts.Colorize))
			colortable = table.getColorTable(opts.Colorize);
		else {
			if (!MatrixIO::parse(opts.Colorize, "color table", colortable))
				return 255;
			assert( matrix.rows() == num_pixels);
		}
	} else {
		colortable = table.getColorTable("blackwhite");
	}
	assert( colortable.cols() == 3 );

//	if (!boost::filesystem::exists(opts.image_file)) {
		// find bounds
		const double min = matrix.minCoeff();
		const double max = matrix.maxCoeff();
		LOG(info, "Data range is [" << min << ":" << max << "]");

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

			 png::rgb_pixel pixel_value =
					 getPixelValue(matrix[i+j*multiplier], min, max, colortable);
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
//	} else {
//		LOG(error, "The given output file exists.");
//		return 1;
//	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	LOG(info, "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start) << ".");

	// exit
	return 0;
}



