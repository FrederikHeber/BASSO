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

	/** Add more sensibility checks.
	 *
	 */
	bool internal_checkSensibility() const;

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
	bool LeftToRight;
	bool BottomToTop;
	bool Colorize;
	bool Flip;
	unsigned int Rotate;

	// secondary values
};

#endif /* MATRIXTOPNGOPTIONS_HPP_ */
