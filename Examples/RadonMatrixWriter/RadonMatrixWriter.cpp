/*
 * RadonMatrixWriter.cpp
 *
 *  Created on: Mar 15, 2016
 *      Author: heber
 */


#include "BassoConfig.h"

#include <Eigen/Dense>

#include <fstream>
#include <iomanip>
#include <iostream>

#include "RadonMatrixWriter/Options/RadonMatrixWriterOptions.hpp"
#include "ComputerTomography/DiscretizedRadon/Backprojection.hpp"

#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"

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

	// set verbosity level
	opts.setVerbosity();

	if (!opts.checkSensibility())
		return 255;
	opts.setSecondaryValues();

	// starting timing
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	BackprojectionMatrix backprojection(
			opts.num_pixel_x,
			opts.num_pixel_y,
			opts.num_angles,
			opts.num_offsets);
	// use transpose of Backprojection as discretized Radon matrox
	Eigen::MatrixXd radon = backprojection.getMatrix().transpose();

	// store matrix to file if desired
	MatrixIO::store(
			opts.radon_matrix.string(),
			"radon matrix",
			radon);

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	BOOST_LOG_TRIVIAL(info) << "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start)
			<< ".";

	// exit
	return 0;
}



