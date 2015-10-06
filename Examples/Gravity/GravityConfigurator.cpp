/*
 * GravityConfigurator.cpp
 *
 *  Created on: Oct 6, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "Log/Logging.hpp"

#include "Options/GravityOptions.hpp"

int main (int argc, char *argv[])
{
	GravityOptions opts;
	opts.init();

	// parse options
	opts.parse(argc, argv);

	if (opts.showHelpConditions(argv[0]))
		return 1;

	// store file to the one given in config
	if (!opts.config_filename.empty()) {
		std::ofstream config_file(opts.config_filename.string().c_str());
		std::cerr << "Writing configuration file "
				<< opts.config_filename.string() << std::endl;
		opts.store(config_file);
	} else {
		std::cerr << "No configuration file specified." << std::endl;
		return 255;
	}

	return 0;
}
