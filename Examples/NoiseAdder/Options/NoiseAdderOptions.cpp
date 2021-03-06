/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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

#include <NoiseAdder/Options/NoiseAdderOptions.hpp>

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
		LOG(debug, "Filename of input set to " << input_file);
	}

	if (vm.count("output")) {
		output_file = vm["output"].as<boost::filesystem::path>();
		LOG(debug, "Writing output vector to " << output_file);
	}

	if (vm.count("noise-level")) {
		noiselevel = vm["noise-level"].as<double>();
		LOG(debug, "Noise level is set to " << noiselevel);
	}

	if (vm.count("relative-level")) {
		relativelevel = vm["relative-level"].as<bool>();
		LOG(debug, "Noise level is set to " << relativelevel);
	}
}

bool NoiseAdderOptions::internal_checkSensibility() const
{
	if (!vm.count("noise-level")) {
		LOG(error, "Noise level is not set");
		return false;
	}

	if ((!vm.count("input")) || (!boost::filesystem::exists(input_file))) {
		LOG(error, "Input file not set or not present.");
		return false;
	}

	if ((!vm.count("output")) || (boost::filesystem::exists(output_file))) {
		LOG(error, "Output file not set or already present.");
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
