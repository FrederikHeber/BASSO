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
			("source-first-factor", po::value< boost::filesystem::path >(),
					"set the file name to parse the first matrix factor")
			("source-second-factor", po::value< boost::filesystem::path >(),
					"set the file name to parse the second matrix factor")
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

	if (vm.count("source-first-factor")) {
		source_first_factor = vm["source-first-factor"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Parsing first matrix factor file from " << source_first_factor.string();
	}

	if (vm.count("source-second-factor")) {
		source_second_factor = vm["source-second-factor"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Parsing second matrix factor file from " << source_second_factor.string();
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
	if (!vm.count("source-first-factor")
			|| !boost::filesystem::exists(source_first_factor)) {
		BOOST_LOG_TRIVIAL(error)
				<< "First matrix factor file not specified or non-existent.";
		return false;
	}
	if (!vm.count("source-second-factor")
			|| !boost::filesystem::exists(source_second_factor)) {
		BOOST_LOG_TRIVIAL(error)
				<< "Second matrix factor file not specified or non-existent.";
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
	writeValue<boost::filesystem::path>(_output, vm,  "source-first-factor");
	writeValue<boost::filesystem::path>(_output, vm,  "source-second-factor");
}
