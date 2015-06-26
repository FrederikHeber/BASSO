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

#include "Options/CommandLineOptions.hpp"

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

	/** Store additional values to output stream.
	 *
	 */
	void internal_store(std::ostream &_output) const;

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
