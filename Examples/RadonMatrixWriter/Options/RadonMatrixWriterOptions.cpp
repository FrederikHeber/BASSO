/*
 * RadonMatrixWriterOptions.cpp
 *
 *  Created on: Mar 15, 2016
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/filesystem.hpp>
#include <iostream>

#include "Log/Logging.hpp"

#include "RadonMatrixWriter/Options/RadonMatrixWriterOptions.hpp"

namespace po = boost::program_options;

RadonMatrixWriterOptions::RadonMatrixWriterOptions() :
		num_pixel_x(0),
		num_pixel_y(0),
		num_angles(0),
		num_offsets(0)
{}

void RadonMatrixWriterOptions::init()
{
	Options::init();

	boost::program_options::options_description desc_radon("RadonMatrixWriter options");

	desc_radon.add_options()
	        ("num-pixels-x", po::value< unsigned int >(),
	        		"set the desired number of pixels in x direction")
			("num-pixels-y", po::value< unsigned int >(),
					"set the desired number of pixels in y direction")
			("num-angles", po::value< unsigned int >(),
					"set the desired number of angle discretization steps")
			("num-offsets", po::value< unsigned int >(),
					"set the desired number of lateral offsets of detector")
			("radon-matrix", po::value< boost::filesystem::path >(),
					"set the file name to write the discretized Radon transformation matrix to")
	        ;

	desc_all.add(desc_radon);
}

void RadonMatrixWriterOptions::parse(int argc, char **argv)
{
	Options::parse(argc,argv);

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

	if (vm.count("num-angles")) {
		num_angles = vm["num-angles"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Number of angle steps was set to " << num_angles;
	}

	if (vm.count("num-offsets")) {
		num_offsets = vm["num-offsets"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Number of offsets steps was set to " << num_offsets;
	}

	if (vm.count("radon-matrix")) {
		radon_matrix = vm["radon-matrix"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of vector was set to " << radon_matrix;
	}
}

bool RadonMatrixWriterOptions::checkSensibility() const
{
	Options::checkSensibility();

	bool status = true;
	status &= checkSensibility_dimensions();

	if (!status)
		showHelpinErrorCase();

	return status;
}

bool RadonMatrixWriterOptions::checkSensibility_dimensions() const
{
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

	if (!vm.count("num-angles")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Number of angle discretization steps not set";
		return false;
	}

	if (!vm.count("num-offsets")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Number of lateral offsets not set";
		return false;
	}

	if ((!vm.count("radon-matrix"))) { // || (boost::filesystem::exists(radon_matrix))) {
		BOOST_LOG_TRIVIAL(error)
				<< "Radon matrix file is not set."; // or already present.";
		return false;
	}

	return true;
}

void RadonMatrixWriterOptions::setSecondaryValues()
{
	Options::setSecondaryValues();
}

void RadonMatrixWriterOptions::store(std::ostream &_output) const
{
	Options::store(_output);

	_output << "# [RadonMatrixWriter]" << std::endl;
	writeValue<unsigned int>(_output, vm,  "num-pixels-x");
	writeValue<unsigned int>(_output, vm,  "num-pixels-y");
	writeValue<unsigned int>(_output, vm,  "num-angles");
	writeValue<unsigned int>(_output, vm,  "num-offsets");
	writeValue<boost::filesystem::path>(_output, vm,  "radon-matrix");
}
