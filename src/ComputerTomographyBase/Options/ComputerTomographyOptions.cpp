/*
 * ComputerTomographyOptions.cpp
 *
 *  Created on: Jul 06, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/filesystem.hpp>
#include <ComputerTomographyBase/Options/ComputerTomographyOptions.hpp>
#include <iostream>

#include "Log/Logging.hpp"

namespace po = boost::program_options;

ComputerTomographyOptions::ComputerTomographyOptions() :
		num_pixel_x(0),
		num_pixel_y(0),
		num_angles(0),
		num_offsets(0)
{}

void ComputerTomographyOptions::internal_init()
{
	boost::program_options::options_description desc_basso("ComputerTomography options");

	desc_basso.add_options()
			("compare-against", po::value< boost::filesystem::path >(),
					"set the file name of the solution to compare against (BregmanDistance)")
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
	        ("rhs", po::value< boost::filesystem::path >(),
	        		"set the vector file of the right-hand side")
			("solution", po::value< boost::filesystem::path >(),
					"set the file name to write solution vector to, i.e. x in y = A*x")
			("solution-image", po::value< boost::filesystem::path >(),
					"set the file name to write image of solution vector to, i.e. A*x")
	        ;

	desc_all.add(desc_basso);
}

void ComputerTomographyOptions::internal_parse()
{
	if (vm.count("compare-against")) {
		comparison_file = vm["compare-against"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Parsing true solution vector from " << comparison_file;
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

	if (vm.count("rhs")) {
		rhs_file = vm["rhs"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of vector was set to " << rhs_file;
	}

	if (vm.count("solution")) {
		solution_file = vm["solution"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Writing solution vector to " << solution_file;
	}

	if (vm.count("solution-image")) {
		solution_image_file = vm["solution-image"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Writing image of solution vector to " << solution_image_file;
	}
}

bool ComputerTomographyOptions::internal_help_conditions() const
{
	if (!vm.count("normx")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Norm of space X normx not set";
		return true;
	}

	if (!vm.count("normy")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Norm of space Y normy not set";
		return true;
	}

	if (!vm.count("num-pixels-x")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Number of pixels in x direction not set";
		return true;
	}

	if (!vm.count("num-pixels-y")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Number of pixels in y direction not set";
		return true;
	}

	if (!vm.count("num-angles")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Number of angle discretization steps not set";
		return true;
	}

	if (!vm.count("num-offsets")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Number of lateral offsets not set";
		return true;
	}

	if ((!vm.count("rhs")) || (!boost::filesystem::exists(rhs_file))) {
		BOOST_LOG_TRIVIAL(error)
				<< "Right-hand side file not set or not present.";
		return true;

	}

	return false;
}

bool ComputerTomographyOptions::internal_checkSensibility() const
{
	return true;

}

void ComputerTomographyOptions::internal_setSecondaryValues()
{}

void ComputerTomographyOptions::internal_store(std::ostream &_output) const
{
	_output << "# [ComputerTomography]" << std::endl;
	writeValue<boost::filesystem::path>(_output, vm,  "compare-against");
	writeValue<unsigned int>(_output, vm,  "num-pixels-x");
	writeValue<unsigned int>(_output, vm,  "num-pixels-y");
	writeValue<unsigned int>(_output, vm,  "num-angles");
	writeValue<unsigned int>(_output, vm,  "num-offsets");
	writeValue<boost::filesystem::path>(_output, vm,  "radon-matrix");
	writeValue<boost::filesystem::path>(_output, vm,  "rhs");
	writeValue<boost::filesystem::path>(_output, vm,  "solution");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-image");
}
