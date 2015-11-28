/*
 * \file NonnegativeMatrixFactors.cpp
 *
 * This will follow [Boutsidis, Gallopoulos, '08].
 *
 *  Created on: Nov 28, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <cmath>
#include <fstream>

#include <boost/chrono.hpp>
#include <Eigen/Dense>

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "Options/NonnegativeMatrixFactorsOptions.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/SpaceElementIO.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"
#include "Solvers/InverseProblemSolver.hpp"


int main (int argc, char *argv[])
{
	// show program information
	showVersion(std::string(argv[0]));
	showCopyright();

	NonnegativeMatrixFactorsOptions opts;
	opts.init();

	// parse options
	opts.parse(argc, argv);

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

	/// parse in source matrix factors
	using namespace MatrixIO;

	Eigen::MatrixXd spectralmatrix;
	if (!MatrixIO::parse(opts.source_first_factor.string(), "first matrix factor", spectralmatrix))
		return 255;
	Eigen::MatrixXd pixelmatrix;
	if (!MatrixIO::parse(opts.source_second_factor.string(), "second matrix factor", pixelmatrix))
		return 255;

	// check dimensionality
	if (spectralmatrix.cols() != pixelmatrix.rows()) {
		BOOST_LOG_TRIVIAL(error)
				<< "Matrix dimensions of both factors do not match";
		return 255;
	}

	/// apply method
	Eigen::MatrixXd nonnegative_spectralmatrix = spectralmatrix;
	Eigen::MatrixXd nonnegative_pixelmatrix = pixelmatrix;

	/// store destination matrix factors
	if (!MatrixIO::store(
			opts.destination_first_factor.string(),
			"first matrix factor",
			nonnegative_spectralmatrix))
		return 255;
	if (!MatrixIO::store(
			opts.destination_second_factor.string(),
			"second matrix factor",
			nonnegative_pixelmatrix))
		return 255;

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	BOOST_LOG_TRIVIAL(info) << "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start)
			<< ".";

	// exit
	return 0;
}

