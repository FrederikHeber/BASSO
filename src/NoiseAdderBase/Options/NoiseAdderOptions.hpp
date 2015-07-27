/*
 * NoiseAdderOptions.hpp
 *
 *  Created on: Jul 27, 2015
 *      Author: heber
 */

#ifndef RANGEPROJECTOROPTIONS_HPP_
#define RANGEPROJECTOROPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>

#include "Options/Options.hpp"

class NoiseAdderOptions : public Options
{
public:
	NoiseAdderOptions();

	/** Add more options to \a desc.
	 *
	 */
	void init();

	/** Adds parsing for additional options.
	 *
	 */
	void parse(int argc, char **argv);

	/** Add more sensibility checks.
	 *
	 */
	bool checkSensibility() const;

	/** Add more secondary values.
	 *
	 */
	void setSecondaryValues();

	/** Store additional values to output stream.
	 *
	 */
	void store(std::ostream &_output) const;

public:

	// primary values
	boost::filesystem::path input_file;
	boost::filesystem::path output_file;
	double noiselevel;
	bool relativelevel;

	// secondary values
};

#endif /* RANGEPROJECTOROPTIONS_HPP_ */
