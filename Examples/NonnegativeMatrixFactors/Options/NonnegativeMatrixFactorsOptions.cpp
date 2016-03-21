/*
 * NonnegativeMatrixFactorsOptions.cpp
 *
 *  Created on: Nov 28, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/filesystem.hpp>
#include <NonnegativeMatrixFactors/Options/NonnegativeMatrixFactorsOptions.hpp>
#include <iostream>

#include "Log/Logging.hpp"

namespace po = boost::program_options;

NonnegativeMatrixFactorsOptions::NonnegativeMatrixFactorsOptions()
{}

void NonnegativeMatrixFactorsOptions::init()
{
	Options::init();

	boost::program_options::options_description desc_nonnegative("NonnegativeMatrixFactors options");

	desc_nonnegative.add_options()
			("matrix", po::value< boost::filesystem::path >(),
					"set the file name to parse the matrix")
			("truncation-dimension", po::value< unsigned int >(),
					"set the truncation dimension")
			("destination-first-factor", po::value< boost::filesystem::path >(),
					"set the file name to write the first nonnegative matrix factor")
			("destination-second-factor", po::value< boost::filesystem::path >(),
					"set the file name to write the second nonnegative matrix factor")
	        ;

	desc_all.add(desc_nonnegative);
}

void NonnegativeMatrixFactorsOptions::parse(int argc, char **argv)
{
	Options::parse(argc,argv);

	if (vm.count("matrix")) {
		matrix = vm["matrix"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Parsing matrix file from " << matrix.string();
	}

	if (vm.count("truncation-dimension")) {
		truncation_dimension = vm["truncation-dimension"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Truncation dimension is " << truncation_dimension;
	}

	if (vm.count("destination-first-factor")) {
		destination_first_factor = vm["destination-first-factor"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Writing first matrix factor file to " << destination_first_factor.string();
	}

	if (vm.count("destination-second-factor")) {
		destination_second_factor = vm["destination-second-factor"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Writing second matrix factor file to " << destination_second_factor.string();
	}
}

bool NonnegativeMatrixFactorsOptions::internal_checkSensibility() const
{
	if (!vm.count("matrix")
			|| !boost::filesystem::exists(matrix)) {
		BOOST_LOG_TRIVIAL(error)
				<< "Matrix file not specified or non-existent.";
		return false;
	}
	if (!vm.count("truncation-dimension")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Truncation_dimension not specified.";
		return false;
	}

	return true;
}

bool NonnegativeMatrixFactorsOptions::checkSensibility() const
{
	bool status = Options::checkSensibility();

	status &= internal_checkSensibility();

	if (!status)
		showHelpinErrorCase();

	return status;
}

void NonnegativeMatrixFactorsOptions::setSecondaryValues()
{
	Options::setSecondaryValues();
}


void NonnegativeMatrixFactorsOptions::store(std::ostream &_output) const
{
	_output << "# [NonnegativeMatrixFactors]" << std::endl;
	writeValue<boost::filesystem::path>(_output, vm,  "destination-first-factor");
	writeValue<boost::filesystem::path>(_output, vm,  "destination-second-factor");
	writeValue<boost::filesystem::path>(_output, vm,  "matrix");
	writeValue<unsigned int>(_output, vm,  "truncation_dimension");
}
