/*
 * MatrixToPNGOptions.cpp
 *
 *  Created on: Oct 09, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include "MatrixToPNG/Colortables/ColorTable.hpp"
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
		Flip(false),
		Rotate(0)
{}

void MatrixToPNGOptions::init()
{
	Options::init();

	boost::program_options::options_description desc_matrixtopng("MatrixToPNG options");

	desc_matrixtopng.add_options()
			("bottom-to-top", po::value< bool >(),
					"set whether to run the columns from bottom to top or reverse")
			("colorize", po::value< std::string >(),
					"use a color table such as redgreen, redblue, bluegreenred or blackwhite (default)")
			("flip", po::value< bool >(),
					"set whether to exchange rows and columns")
			("image", po::value< boost::filesystem::path >(),
					"set the file name of output image")
			("left-to-right", po::value< bool >(),
					"set whether to run the rows from left to right or reverse")
			("matrix", po::value< boost::filesystem::path >(),
					"set the file name of input matrix")
	        ("num-pixels-x", po::value< unsigned int >(),
	        		"set the desired number of pixels in x direction")
			("num-pixels-y", po::value< unsigned int >(),
					"set the desired number of pixels in y direction")
			("rotate", po::value< unsigned int >(),
					"set final rotation to none(0), -90 (1), -180(2), +90 (3)")
	        ;

	desc_all.add(desc_matrixtopng);
}

void MatrixToPNGOptions::parse(int argc, char **argv)
{
	Options::parse(argc,argv);

	if (vm.count("bottom-to-top")) {
		BottomToTop = vm["bottom-to-top"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "We go through columns from "
			<< (BottomToTop ? "bottom to top" : "top to bottom");
	}

	if (vm.count("colorize")) {
		Colorize = vm["colorize"].as<std::string>();
		if (!Colorize.empty())
			BOOST_LOG_TRIVIAL(debug)
				<< "We do use colors to designate positive and negative areas, using "
				<< Colorize << ".";
		else
			BOOST_LOG_TRIVIAL(debug)
				<< "We don't use colors to designate positive and negative areas.";
	}

	if (vm.count("flip")) {
		Flip = vm["flip"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "We " << (Flip ? "do" : "don't") << " exchange rows and columns";
	}

	if (vm.count("image")) {
		image_file = vm["image"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Writing image from " << image_file;
	}

	if (vm.count("left-to-right")) {
		LeftToRight = vm["left-to-right"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "We go through rows from "
			<< (LeftToRight ? "left to right" : "right to left");
	}

	if (vm.count("matrix")) {
		matrix_file = vm["matrix"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Parsing matrix from " << matrix_file;
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

	if (vm.count("rotate")) {
		Rotate = vm["rotate"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "We rotate to " << Rotate;
	}
}

bool MatrixToPNGOptions::internal_checkSensibility() const
{
	if (!vm.count("matrix")) {
		BOOST_LOG_TRIVIAL(error)
				<< "matrix filename not set";
		return false;
	}

	if (vm.count("colorize")) {
		// check whether it's a file name or a key in the color table
		ColorTable table;
		boost::filesystem::path filepath(Colorize);
		if ((!table.isKeyPresent(Colorize)) && (!boost::filesystem::exists(filepath))) {
			BOOST_LOG_TRIVIAL(error)
				<< "Either the color table name is mis-spelled or the given color matrix filename does not exist";
			return false;
		}
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

	if (vm.count("rotate"))
		status &= (Rotate >= 0) && (Rotate <=3);

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
	writeValue<bool>(_output, vm,  "bottom-to-top");
	writeValue<std::string>(_output, vm,  "colorize");
	writeValue<bool>(_output, vm,  "flip");
	writeValue<boost::filesystem::path>(_output, vm,  "image");
	writeValue<bool>(_output, vm,  "left-to-right");
	writeValue<boost::filesystem::path>(_output, vm,  "matrix");
	writeValue<unsigned int>(_output, vm,  "num-pixels-x");
	writeValue<unsigned int>(_output, vm,  "num-pixels-y");
	writeValue<unsigned int>(_output, vm,  "rotate");
}
