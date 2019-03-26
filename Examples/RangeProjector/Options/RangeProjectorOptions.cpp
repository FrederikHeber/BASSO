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
 * RangeProjectorOptions.cpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/filesystem.hpp>
#include <RangeProjector/Options/RangeProjectorOptions.hpp>
#include <iostream>

#include "Log/Logging.hpp"

namespace po = boost::program_options;

void RangeProjectorOptions::internal_init()
{
	boost::program_options::options_description desc_basso("RangeProjector options");

	desc_basso.add_options()
	        ("matrix", po::value< boost::filesystem::path >(),
	        		"set the forward operator matrix file")
	        ("rhs", po::value< boost::filesystem::path >(),
	        		"set the vector file of the right-hand side")
			("solution", po::value< boost::filesystem::path >(),
					"set the file name to write solution vector to, i.e. x in y = A*x")
			("solution-image", po::value< boost::filesystem::path >(),
					"set the file name to write image of solution vector to, i.e. A*x")
	        ;

	desc_all.add(desc_basso);
}

void RangeProjectorOptions::internal_parse()
{
	if (vm.count("matrix")) {
		matrix_file = vm["matrix"].as<boost::filesystem::path>();
		LOG(debug, "Filename of matrix was set to " << matrix_file);
	}

	if (vm.count("rhs")) {
		rhs_file = vm["rhs"].as<boost::filesystem::path>();
		LOG(debug, "Filename of vector was set to " << rhs_file);
	}

	if (vm.count("solution")) {
		solution_file = vm["solution"].as<boost::filesystem::path>();
		LOG(debug, "Writing solution vector to " << solution_file);
	}

	if (vm.count("solution-image")) {
		solution_image_file = vm["solution-image"].as<boost::filesystem::path>();
		LOG(debug, "Writing image of solution vector to " << solution_image_file);
	}
}

bool RangeProjectorOptions::internal_checkSensibility() const
{
	if ((!vm.count("matrix")) || (!boost::filesystem::exists(matrix_file))) {
		LOG(error, "Matrix file not set or not present.");
		return false;

	}

	if ((!vm.count("rhs")) || (!boost::filesystem::exists(rhs_file))) {
		LOG(error, "Right-hand side file not set or not present.");
		return false;

	}

	return true;
}

void RangeProjectorOptions::internal_setSecondaryValues()
{}

void RangeProjectorOptions::internal_store(std::ostream &_output) const
{
	_output << "# [RangeProjector]" << std::endl;
	writeValue<boost::filesystem::path>(_output, vm,  "matrix");
	writeValue<boost::filesystem::path>(_output, vm,  "rhs");
	writeValue<boost::filesystem::path>(_output, vm,  "solution");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-image");
}
