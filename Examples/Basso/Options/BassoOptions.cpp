/*
 * BassoOptions.cpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/filesystem.hpp>
#include <Basso/Options/BassoOptions.hpp>

#include <iostream>
#include <sstream>

#include "Log/Logging.hpp"

namespace po = boost::program_options;

void BassoOptions::internal_init()
{
	boost::program_options::options_description desc_basso("Basso options");

	desc_basso.add_options()
			("compare-against", po::value< boost::filesystem::path >(),
					"set the file name of the solution to compare against (BregmanDistance)")
	        ("matrix", po::value< std::vector<boost::filesystem::path> >()->multitoken(),
	        		"set the forward operator matrix file, matrix factors if multiple files are given")
	        ("rhs", po::value< boost::filesystem::path >(),
	        		"set the vector file of the right-hand side")
			("solution", po::value< boost::filesystem::path >(),
					"set the file name to write solution vector to, i.e. x in y = A*x")
			("solution-image", po::value< boost::filesystem::path >(),
					"set the file name to write image of solution vector to, i.e. A*x")
	        ;

	desc_all.add(desc_basso);
}

void BassoOptions::internal_parse()
{
	if (vm.count("compare-against")) {
		comparison_file = vm["compare-against"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Parsing true solution vector from " << comparison_file;
	}

	if (vm.count("matrix")) {
		matrix_file = vm["matrix"].as< std::vector<boost::filesystem::path> >();
		std::stringstream output;
		output << matrix_file;
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of matrix was set to " << output.str();
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

bool BassoOptions::internal_checkSensibility() const
{
	// check whether both or none are given
	switch (matrix_file.size())
	{
	case 1:
		if (!boost::filesystem::exists(matrix_file[0])) {
			BOOST_LOG_TRIVIAL(error)
					<< "File of matrix does not exist.";
			return false;
		}
		break;
	case 2:
		if ((!boost::filesystem::exists(matrix_file[0]))
				|| (!boost::filesystem::exists(matrix_file[1]))) {
			BOOST_LOG_TRIVIAL(error)
					<< "At least one of the matrix files does not exist.";
			return false;
		}
		break;
	default:
		BOOST_LOG_TRIVIAL(error)
				<< "More than two matrix files given.";
		return false;
		break;
	}

	if ((!vm.count("rhs")) || (!boost::filesystem::exists(rhs_file))) {
		BOOST_LOG_TRIVIAL(error)
				<< "Right-hand side file not set or not present.";
		return false;

	}

	return true;
}

void BassoOptions::internal_setSecondaryValues()
{}

void BassoOptions::internal_store(std::ostream &_output) const
{
	_output << "# [Basso]" << std::endl;
	writeValue<boost::filesystem::path>(_output, vm,  "compare-against");
	if (matrix_file.size() != 0) {
		for (std::vector<boost::filesystem::path>::const_iterator iter = matrix_file.begin();
				iter != matrix_file.end();++iter)
			_output << "\tmatrix = " << *iter << std::endl;
	}
	writeValue<boost::filesystem::path>(_output, vm,  "rhs");
	writeValue<boost::filesystem::path>(_output, vm,  "solution");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-image");
}
