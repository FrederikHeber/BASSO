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
 * NoiseAdder.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <string>

#include <boost/assign.hpp>

#include "NoiseAdder/Options/NoiseAdderOptions.hpp"

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"

using namespace boost::assign;

int main (int argc, char *argv[])
{
	// show program information
	showVersion(std::string(argv[0]));
	showCopyright();

	NoiseAdderOptions opts;
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

	// parse matrix and vector files into instances
	Eigen::MatrixXd input_element;
	{
		using namespace MatrixIO;

		{
			std::ifstream ist(opts.input_file.string().c_str());
			if (ist.good())
				try {
					ist >> input_element;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse rhs from " << opts.input_file.string() << std::endl;
					return 255;
				}
			else {
				std::cerr << "Failed to open " << opts.input_file.string() << std::endl;
				return 255;
			}
		}
	}

	// add random noise to matrix
	Eigen::MatrixXd output_element(input_element.rows(), input_element.cols());
	output_element.setRandom();
	if (opts.relativelevel) {
		const double maximum = input_element.maxCoeff();
		output_element *= maximum*opts.noiselevel;
	} else
		output_element *= opts.noiselevel;
	output_element += input_element;

	// writing solution
	{
		using namespace MatrixIO;
		if (!opts.output_file.string().empty()) {
			std::ofstream ost(opts.output_file.string().c_str());
			if (ost.good())
				try {
					ost << std::setprecision(10) << output_element;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write solution to file.\n";
				}
			else {
				std::cerr << "Failed to open " << opts.output_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No solution file given." << std::endl;
		}
	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	LOG(info, "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start) << ".");

	// exit
	return 0;
}



