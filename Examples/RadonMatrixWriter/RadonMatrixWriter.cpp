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
	LOG(info, "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start) << ".");

	// exit
	return 0;
}



