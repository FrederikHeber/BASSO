/*
 * MatrixToPNGOptions.hpp
 *
 *  Created on: Oct 09, 2015
 *      Author: heber
 */

#ifndef MATRIXTOPNGOPTIONS_HPP_
#define MATRIXTOPNGOPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>

#include "Options/Options.hpp"

class MatrixToPNGOptions : public Options
{
public:
	MatrixToPNGOptions();

public:
	/** Add more options to \a desc.
	 *
	 */
	void init();

	/** Adds parsing for additional options.
	 *
	 */
	void parse(int argc, char **argv);

	/** Add more conditions when help is shown.
	 *
	 */
	bool help_conditions() const;

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
	boost::filesystem::path matrix_file;
	boost::filesystem::path image_file;
	unsigned int num_pixel_x;
	unsigned int num_pixel_y;

	// secondary values
};

#endif /* MATRIXTOPNGOPTIONS_HPP_ */
