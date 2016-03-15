/*
 * RadonMatrixWriterConfigurator.cpp
 *
 *  Created on: Mar 15, 2016
 *      Author: heber
 */


#include "BassoConfig.h"

#include <boost/exception/all.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <fstream>

#include "Log/Logging.hpp"

#include "RadonMatrixWriter/Options/RadonMatrixWriterOptions.hpp"

int main (int argc, char *argv[])
{
	// show program information
	showVersion(std::string(argv[0]));
	showCopyright();

	RadonMatrixWriterOptions opts;
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
