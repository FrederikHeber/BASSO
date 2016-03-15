/*
 * \file NonnegativeMatrixFactors.cpp
 *
 * This follows [Boutsidis, Gallopoulos, '08].
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

static Eigen::VectorXd getPositiveSection(const Eigen::VectorXd &_vector)
{
	Eigen::VectorXd returnvector = _vector.array().abs();
	returnvector += _vector;
	returnvector *= .5;
	return returnvector;
}

static Eigen::VectorXd getNegativeSection(const Eigen::VectorXd &_vector)
{
	Eigen::VectorXd returnvector = (-1.*_vector).array().abs();
	returnvector -= _vector;
	returnvector *= -.5;
	return returnvector;
}

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
	const size_t inner_dimension = spectralmatrix.cols();

	/// apply method
	// create destination matrices
	Eigen::MatrixXd nonnegative_spectralmatrix =
			Eigen::MatrixXd(spectralmatrix.rows(), spectralmatrix.cols());
	Eigen::MatrixXd nonnegative_pixelmatrix =
			Eigen::MatrixXd(pixelmatrix.rows(), pixelmatrix.cols());

	// loop through other columns/row and extract positive factor in rank-2 decomp.
	for (size_t index = 0; index < inner_dimension; ++index) {
		Eigen::VectorXd x = spectralmatrix.col(index);
		Eigen::VectorXd y = pixelmatrix.row(index);
		Eigen::VectorXd xp = getPositiveSection(x);
		Eigen::VectorXd xn = getNegativeSection(x);
		Eigen::VectorXd yp = getPositiveSection(y);
		Eigen::VectorXd yn = getNegativeSection(y);
		const double xp_norm = xp.stableNorm();
		const double xn_norm = xn.stableNorm();
		const double yp_norm = yp.stableNorm();
		const double yn_norm = yn.stableNorm();
		const double mp = xp_norm * yp_norm;
		const double mn = xn_norm * yn_norm;
		assert( mp >= 0. );
		assert( mn >= 0. );
		if ((fabs(mp) > BASSOTOLERANCE) || (fabs(mn) > BASSOTOLERANCE)) {
			if (mp > mn) {
				nonnegative_spectralmatrix.col(index) = (mp/xp_norm)*xp;
				nonnegative_pixelmatrix.row(index) = (mp/yp_norm)*yp;
			} else {
				nonnegative_spectralmatrix.col(index) = (-mn/xn_norm)*xn;
				nonnegative_pixelmatrix.row(index) = (-mn/yn_norm)*yn;
			}
		} else {
			nonnegative_spectralmatrix.col(index).setZero();
			nonnegative_pixelmatrix.row(index).setZero();
		}
	}

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

