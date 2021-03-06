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
 * Options.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Options.hpp"

#include <fstream>

#include "Log/Logging.hpp"

Options::Options() :
	desc_all("Configuration Options"),
	verbose(0)
{}

namespace po = boost::program_options;

void Options::init()
{
	boost::program_options::options_description desc_run("Run options");

	desc_config.add_options()
			("config", po::value<boost::filesystem::path>(),
					"filename of a configuration file containing default option parameters. Note that other command-line values override those in this file.")
			;

	desc_run.add_options()
			("help",
					"produce help message")
			("verbose", po::value<unsigned int>(),
					"set the amount of verbosity")
			("version",
					"show version information")
					;

	desc_all
		.add(desc_config)
		.add(desc_run);

	setVerbosity();
}

void Options::parse(int argc, char **argv)
{
	/// parse the usual command line options
	po::store(po::parse_command_line(argc, argv, desc_all), vm);
	po::notify(vm);

	/// additionally parse config file
	if (vm.count("config")) {
		config_filename = vm["config"].as<boost::filesystem::path>();
		std::ifstream config_file(config_filename.string().c_str());
		po::store(po::parse_config_file(config_file, desc_all), vm);
		po::notify(vm);
	}

	if (vm.count("verbose")) {
		verbose = vm["verbose"].as<unsigned int>();
		LOG(debug, "Verbose set to " << verbose);
		setVerbosity();
	}
}

bool Options::showHelpConditions(const char * const program_name) const
{
	if (vm.count("version"))
		return true;
	if (vm.count("help")) {
	    std::cout << desc_all << "\n";
	    return true;
	} else {
		return false;
	}
}

void Options::showHelpinErrorCase() const
{
	std::cout << "There was an error parsing options, use '--help' to learn more." << std::endl;
}

void Options::setVerbosity() const
{
	Verbose::setVerbosity(verbose);
}

void Options::store(std::ostream &_output) const
{
	_output << "# [Run]" << std::endl;
	writeValue<unsigned int>(_output, vm,  "verbose");
}
