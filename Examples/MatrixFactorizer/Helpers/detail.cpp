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
	_opts.parse(_argc, _argv);

	if (_opts.showHelpConditions(_argv[0]))
		return 1;

	// set verbosity level
	_opts.setVerbosity();

	if (!_opts.checkSensibility())
		return 255;
	_opts.setSecondaryValues();

	return 0;
}

int detail::parseDataFile(
		const std::string &_filename,
		Eigen::MatrixXd &_data)
{
	using namespace MatrixIO;

	if (!MatrixIO::parse(_filename, "data matrix", _data))
		return 255;

	// print parsed matrix and vector if small or high verbosity requested
	if ((_data.innerSize() > 10) || (_data.outerSize() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
			<< "We solve for Y=K*X with Y =\n"
			<< _data << "." << std::endl;
	} else {
		BOOST_LOG_TRIVIAL(info)
					<< "We solve for Y=K*X with Y =\n"
					<< _data << "." << std::endl;
	}

	return 0;
}

void detail::constructStartingMatrices(
		Eigen::MatrixXd &_spectral_matrix,
		Eigen::MatrixXd &_pixel_matrix,
		const bool _nonnegative
		)
{
	TraceRenormalizer renormalizer;
	// generate random matrix drawn from [-1,1]
	_spectral_matrix.setRandom();
	// push to [0,1] if non-negativity desired
	if (_nonnegative) {
		_spectral_matrix += Eigen::MatrixXd::Ones(_spectral_matrix.rows(), _spectral_matrix.cols());
		_spectral_matrix *= .5;
	}
	renormalizer.renormalize(_spectral_matrix);
	if ((_spectral_matrix.innerSize() > 10) || (_spectral_matrix.outerSize() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
				<< "Initial spectral matrix is\n" << _spectral_matrix;
	} else {
		BOOST_LOG_TRIVIAL(info)
				<< "Initial spectral matrix is\n" << _spectral_matrix;
	}
	// set second matrix to zero
	_pixel_matrix.setZero();
}

