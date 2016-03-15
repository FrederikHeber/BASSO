/*
 * RadonMatrixWriterOptions.hpp
 *
 *  Created on: Mar 15, 2016
 *      Author: heber
 */

#ifndef RADONMATRIXWRITEROPTIONS_HPP_
#define RADONMATRIXWRITEROPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>

#include "Options/Options.hpp"

class RadonMatrixWriterOptions : public Options
{
public:
	RadonMatrixWriterOptions();

	/** Initializes the options.
	 *
	 */
	void init();

	/** Performs the actual parsing.
	 *
	 * \param argc argument count
	 * \param **argv array of arguments
	 */
	void parse(int argc, char **argv);

	/** Checks whether the values make any sense.
	 *
	 * \return true - values sensible, false - not
	 */
	bool checkSensibility() const;

	/** Sets other values that require checkSensibility().
	 *
	 */
	void setSecondaryValues();

	/** Stores all variables as "key = value" pairs per line under the given
	 * stream \a _output.
	 *
	 * @param _output output stream to write to
	 */
	void store(std::ostream &_output) const;

private:
	bool checkSensibility_dimensions() const;

public:

	// primary values
	unsigned int num_pixel_x;
	unsigned int num_pixel_y;
	unsigned int num_angles;
	unsigned int num_offsets;
	boost::filesystem::path radon_matrix;

	// secondary values
};

#endif /* RADONMATRIXWRITEROPTIONS_HPP_ */
