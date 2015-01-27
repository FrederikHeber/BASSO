/*
 * BassoOptions.hpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#ifndef BASSOOPTIONS_HPP_
#define BASSOOPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>

#include "CommandLineOptions/CommandLineOptions.hpp"

class BassoOptions : public CommandLineOptions
{
public:
	BassoOptions();

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

	/** Add more secondary values.
	 *
	 */
	void internal_setSecondaryValues();

public:

	// primary values
	boost::filesystem::path comparison_file;
	boost::filesystem::path matrix_file;
	unsigned int maxiter;
	double maxwalltime;
	boost::filesystem::path rhs_file;
	boost::filesystem::path solution_file;
	boost::filesystem::path solution_image_file;

	// secondary values
};

#endif /* BASSOOPTIONS_HPP_ */
