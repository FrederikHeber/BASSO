/*
 * MatrixFactorizerOptions.hpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZEROPTIONS_HPP_
#define MATRIXFACTORIZEROPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>

#include "CommandLineOptions/CommandLineOptions.hpp"

class MatrixFactorizerOptions : public CommandLineOptions
{
public:
	MatrixFactorizerOptions();

private:
	/** Add more options to \a desc.
	 *
	 */
	void internal_init();

	/** Adds parsing for additional options.
	 *
	 */
	void internal_parse();

	/** Add more conditions when help is shown.
	 *
	 */
	bool internal_help_conditions() const;

	/** Add more sensibility checks.
	 *
	 */
	bool internal_checkSensibility() const;

public:
	boost::filesystem::path data_file;
	unsigned int inner_iterations;
	unsigned int max_loops;
	boost::filesystem::path solution_factor_one_file;
	boost::filesystem::path solution_factor_two_file;
	unsigned int sparse_dim;
};

#endif /* MATRIXFACTORIZEROPTIONS_HPP_ */
