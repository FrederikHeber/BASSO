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
		BOOST_LOG_TRIVIAL(trace)
			<< _name << " is "
			<< _data << "." << std::endl;
	} else {
		BOOST_LOG_TRIVIAL(info)
					<< _name << " is "
					<< _data << "." << std::endl;
	}

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
		if (_matrix.innerSize() > _matrix.outerSize())
			BOOST_LOG_TRIVIAL(trace)
					<< "Initial transposed matrix is\n" << _matrix.transpose();
		else
			BOOST_LOG_TRIVIAL(trace)
					<< "Initial matrix is\n" << _matrix;
	} else {
		if (_matrix.innerSize() > _matrix.outerSize())
			BOOST_LOG_TRIVIAL(info)
					<< "Initial transposed matrix is\n" << _matrix.transpose();
		else
			BOOST_LOG_TRIVIAL(info)
					<< "Initial matrix is\n" << _matrix;
	}

}

void detail::constructZeroMatrix(
		Eigen::MatrixXd &_matrix
		)
{
	_matrix.setZero();
}

