/*
 * GravityOptions.hpp
 *
 *  Created on: Oct 06, 2015
 *      Author: heber
 */

#ifndef GRAVITYOPTIONS_HPP_
#define GRAVITYOPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>

#include "Options/CommandLineOptions.hpp"

class GravityOptions : public CommandLineOptions
{
public:
	GravityOptions();

private:
	/** Add more options to \a desc.
	 *
	 */
	void internal_init();

	/** Adds parsing for additional options.
	 *
	 */
	void internal_parse();

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
	double depth;
	unsigned int discretization;
	boost::filesystem::path gravityfield_file;
	boost::filesystem::path density_file;

	// secondary values
};

#endif /* GRAVITYOPTIONS_HPP_ */
