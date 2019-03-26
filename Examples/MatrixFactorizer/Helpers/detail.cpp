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
 * detail.cpp
 *
 *  Created on: Jun 29, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "detail.hpp"

#include <fstream>
#include <string>

#include "Log/Logging.hpp"
#include "MatrixFactorizer/Options/MatrixFactorizerOptions.hpp"
#include "MatrixFactorizer/ScalingAmbiguityRemover/TraceRenormalizer.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "MatrixIO/MatrixIOExceptions.hpp"

int detail::parseOptions(
		int _argc,
		char ** _argv,
		MatrixFactorizerOptions &_opts)
{
	_opts.init();

	// parse options
	try {
		_opts.parse(_argc, _argv);
	} catch (std::exception &e) {
		std::cerr << "An error occurred: "
				<< e.what()
				<< std::endl;
		return 255;
	}

	if (_opts.showHelpConditions(_argv[0]))
		return 1;

	// set verbosity level
	_opts.setVerbosity();

	if (!_opts.checkSensibility())
		return 255;
	_opts.setSecondaryValues();

	return 0;
}

int detail::parseFactorFile(
		const std::string &_filename,
		Eigen::MatrixXd &_data,
		const std::string &_name)
{
	using namespace MatrixIO;

	if (!MatrixIO::parse(_filename, _name, _data))
		return 255;

	// print parsed matrix and vector if small or high verbosity requested
	if ((_data.innerSize() > 10) || (_data.outerSize() > 10)) {
		LOG(trace, _name << " is " << _data << "." << std::endl);
	} else {
		LOG(info, _name << " is " << _data << "." << std::endl);
	}

	return 0;
}

void detail::constructRandomMatrix(
		Eigen::MatrixXd &_matrix
		)
{
	TraceRenormalizer renormalizer;
	// generate random matrix drawn from [-1,1]
	_matrix.setRandom();
	// push to [0,1] if non-negativity desired
	renormalizer.renormalize(_matrix);
	if ((_matrix.innerSize() > 10) || (_matrix.outerSize() > 10)) {
		if (_matrix.innerSize() > _matrix.outerSize()) {
			LOG(trace, "Initial transposed matrix is\n" << _matrix.transpose());
		} else {
			LOG(trace, "Initial matrix is\n" << _matrix);
		}
	} else {
		if (_matrix.innerSize() > _matrix.outerSize()) {
			LOG(info, "Initial transposed matrix is\n" << _matrix.transpose());
		} else {
			LOG(info, "Initial matrix is\n" << _matrix);
		}
	}

}

void detail::constructZeroMatrix(
		Eigen::MatrixXd &_matrix
		)
{
	_matrix.setZero();
}

