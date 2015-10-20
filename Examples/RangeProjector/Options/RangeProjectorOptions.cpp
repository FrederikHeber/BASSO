/*
 * RangeProjectorOptions.cpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/filesystem.hpp>
#include <RangeProjector/Options/RangeProjectorOptions.hpp>
#include <iostream>

#include "Log/Logging.hpp"

namespace po = boost::program_options;

void RangeProjectorOptions::internal_init()
{
	boost::program_options::options_description desc_basso("RangeProjector options");

	desc_basso.add_options()
	        ("matrix", po::value< boost::filesystem::path >(),
	        		"set the forward operator matrix file")
	        ("rhs", po::value< boost::filesystem::path >(),
	        		"set the vector file of the right-hand side")
			("solution", po::value< boost::filesystem::path >(),
					"set the file name to write solution vector to, i.e. x in y = A*x")
			("solution-image", po::value< boost::filesystem::path >(),
					"set the file name to write image of solution vector to, i.e. A*x")
	        ;

	desc_all.add(desc_basso);
}

void RangeProjectorOptions::internal_parse()
{
	if (vm.count("matrix")) {
		matrix_file = vm["matrix"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of matrix was set to " << matrix_file;
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

bool RangeProjectorOptions::internal_checkSensibility() const
{
	if ((!vm.count("matrix")) || (!boost::filesystem::exists(matrix_file))) {
		BOOST_LOG_TRIVIAL(error)
				<< "Matrix file not set or not present.";
		return false;

	}

	if ((!vm.count("rhs")) || (!boost::filesystem::exists(rhs_file))) {
		BOOST_LOG_TRIVIAL(error)
				<< "Right-hand side file not set or not present.";
		return false;

	}

	return true;
}

void RangeProjectorOptions::internal_setSecondaryValues()
{}

void RangeProjectorOptions::internal_store(std::ostream &_output) const
{
	_output << "# [RangeProjector]" << std::endl;
	writeValue<boost::filesystem::path>(_output, vm,  "matrix");
	writeValue<boost::filesystem::path>(_output, vm,  "rhs");
	writeValue<boost::filesystem::path>(_output, vm,  "solution");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-image");
}
