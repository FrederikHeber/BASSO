/*
 * MatrixToPNGOptions.cpp
 *
 *  Created on: Oct 09, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include "MatrixToPNG/Options/MatrixToPNGOptions.hpp"

#include <boost/filesystem.hpp>
#include <iostream>

#include "Log/Logging.hpp"

namespace po = boost::program_options;

MatrixToPNGOptions::MatrixToPNGOptions() :
		num_pixel_x(0),
		num_pixel_y(0),
		LeftToRight(false),
		BottomToTop(false),
		Flip(false)
{}

void MatrixToPNGOptions::init()
{
	Options::init();

	boost::program_options::options_description desc_matrixtopng("MatrixToPNG options");

	desc_matrixtopng.add_options()
			("matrix", po::value< boost::filesystem::path >(),
					"set the file name of input matrix")
			("image", po::value< boost::filesystem::path >(),
					"set the file name of output image")
	        ("num-pixels-x", po::value< unsigned int >(),
	        		"set the desired number of pixels in x direction")
			("num-pixels-y", po::value< unsigned int >(),
					"set the desired number of pixels in y direction")
			("left-to-right", po::value< bool >(),
					"set whether to run the rows from left to right or reverse")
			("bottom-to-top", po::value< bool >(),
					"set whether to run the columns from bottom to top or reverse")
			("flip", po::value< bool >(),
					"set whether to exchange rows and columns")
	        ;

	desc_all.add(desc_matrixtopng);
}

void MatrixToPNGOptions::parse(int argc, char **argv)
{
	Options::parse(argc,argv);

	if (vm.count("matrix")) {
		matrix_file = vm["matrix"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Parsing matrix from " << matrix_file;
	}

	if (vm.count("image")) {
		image_file = vm["image"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Writing image from " << image_file;
	}

	if (vm.count("num-pixels-x")) {
		num_pixel_x = vm["num-pixels-x"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Number of x pixels was set to " << num_pixel_x;
	}

	if (vm.count("num-pixels-y")) {
		num_pixel_y = vm["num-pixels-y"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Number of y pixels was set to " << num_pixel_y;
	}

	if (vm.count("left-to-right")) {
		LeftToRight = vm["left-to-right"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "We go through rows from "
			<< (LeftToRight ? "left to right" : "right to left");
	}

	if (vm.count("bottom-to-top")) {
		BottomToTop = vm["bottom-to-top"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "We go through columns from "
			<< (BottomToTop ? "bottom to top" : "top to bottom");
	}


	if (vm.count("flip")) {
		Flip = vm["flip"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "We " << (Flip ? "do" : "don't") << " exchange rows and columns";
	}
}

bool MatrixToPNGOptions::internal_checkSensibility() const
{
	if (!vm.count("matrix")) {
		BOOST_LOG_TRIVIAL(error)
				<< "matrix filename not set";
		return false;
	}

	if (!vm.count("image")) {
		BOOST_LOG_TRIVIAL(error)
				<< "image filename not set";
		return false;
	}

	if (!vm.count("num-pixels-x")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Number of pixels in x direction not set";
		return false;
	}

	if (!vm.count("num-pixels-y")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Number of pixels in y direction not set";
		return false;
	}

	return true;
}

bool MatrixToPNGOptions::checkSensibility() const
{
	bool status = Options::checkSensibility();

	status &= internal_checkSensibility();

	if (!status)
		showHelpinErrorCase();

	return status;
}

void MatrixToPNGOptions::setSecondaryValues()
{
	Options::setSecondaryValues();
}

void MatrixToPNGOptions::store(std::ostream &_output) const
{
	Options::store(_output);

	_output << "# [MatrixToPNG]" << std::endl;
	writeValue<boost::filesystem::path>(_output, vm,  "matrix");
	writeValue<boost::filesystem::path>(_output, vm,  "image");
	writeValue<unsigned int>(_output, vm,  "num-pixels-x");
	writeValue<unsigned int>(_output, vm,  "num-pixels-y");
	writeValue<bool>(_output, vm,  "left-to-right");
	writeValue<bool>(_output, vm,  "bottom-to-top");
	writeValue<bool>(_output, vm,  "flip");
}
