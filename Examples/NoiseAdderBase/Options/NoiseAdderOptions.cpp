/*
 * NoiseAdderOptions.cpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include <iostream>

#include <boost/filesystem.hpp>

#include "Log/Logging.hpp"

#include <NoiseAdderBase/Options/NoiseAdderOptions.hpp>

namespace po = boost::program_options;

NoiseAdderOptions::NoiseAdderOptions() :
		noiselevel(0.1),
		relativelevel(false)
{}

void NoiseAdderOptions::init()
{
	Options::init();

	boost::program_options::options_description desc_run("NoiseAdder options");

	desc_run.add_options()
	        ("input", po::value< boost::filesystem::path >(),
	        		"set the file name of the vector/matrix to be perturbed")
			("output", po::value< boost::filesystem::path >(),
					"set the file name to write pertured output to")
			("noise-level", po::value< double >(),
					"set noise level to perturb given input with")
			("relative-level", po::value< bool >(),
					"set whether to add the noise as absolute (false) or relative to maximum component (true) with respect to noise level")
	        ;

	desc_all.add(desc_run);
}

void NoiseAdderOptions::parse(int argc, char **argv)
{
	Options::parse(argc,argv);

	if (vm.count("input")) {
		input_file = vm["input"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of input set to " << input_file;
	}

	if (vm.count("output")) {
		output_file = vm["output"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Writing output vector to " << output_file;
	}

	if (vm.count("noise-level")) {
		noiselevel = vm["noise-level"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Noise level is set to " << noiselevel;
	}

	if (vm.count("relative-level")) {
		relativelevel = vm["relative-level"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Noise level is set to " << relativelevel;
	}
}

bool NoiseAdderOptions::internal_checkSensibility() const
{
	if (!vm.count("noise-level")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Noise level is not set";
		return false;
	}

	if ((!vm.count("input")) || (!boost::filesystem::exists(input_file))) {
		BOOST_LOG_TRIVIAL(error)
				<< "Input file not set or not present.";
		return false;
	}

	if ((!vm.count("output")) || (boost::filesystem::exists(output_file))) {
		BOOST_LOG_TRIVIAL(error)
				<< "Output file not set or already present.";
		return false;
	}

	return true;
}

bool NoiseAdderOptions::checkSensibility() const
{
	bool status = Options::checkSensibility();

	status &= internal_checkSensibility();

	if (!status)
		showHelpinErrorCase();

	return status;
}

void NoiseAdderOptions::setSecondaryValues()
{
	Options::setSecondaryValues();
}

void NoiseAdderOptions::store(std::ostream &_output) const
{
	Options::store(_output);

	_output << "# [NoiseAdder]" << std::endl;
	writeValue<boost::filesystem::path>(_output, vm,  "input");
	writeValue<boost::filesystem::path>(_output, vm,  "output");
	writeValue<double>(_output, vm,  "noise-level");
	writeValue<bool>(_output, vm,  "relative-level");
}
