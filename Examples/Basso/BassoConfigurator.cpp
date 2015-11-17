/*
 * BassoConfigurator.cpp
 *
 *  Created on: Apr 7, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "Log/Logging.hpp"

#include <Basso/Options/BassoOptions.hpp>

int main (int argc, char *argv[])
{
	BassoOptions opts;
	opts.init();

	// parse options
	try {
		opts.parse(argc, argv);
	} catch (std::exception &e) {
		std::cerr << "An error occurred: "
				<< e.what()
				<< std::endl;
		return 255;
	}


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
