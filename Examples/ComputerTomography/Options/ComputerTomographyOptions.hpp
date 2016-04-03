/*
 * ComputerTomographyOptions.hpp
 *
 *  Created on: Jul 06, 2015
 *      Author: heber
 */

#ifndef COMPUTERTOMOGRAPHYOPTIONS_HPP_
#define COMPUTERTOMOGRAPHYOPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>

#include <vector>

#include "Options/CommandLineOptions.hpp"

class ComputerTomographyOptions : public CommandLineOptions
{
public:
	ComputerTomographyOptions();

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

private:

	/** Helper function to check validity of radon matrix file parameters.
	 *
	 */
	bool internal_checkSensibility_radonfactors() const;

public:

	// primary values
	boost::filesystem::path comparison_file;
	double noiselevel;
	boost::filesystem::path noisy_sinogram_file;
	unsigned int num_pixel_x;
	unsigned int num_pixel_y;
	unsigned int num_angles;
	unsigned int num_offsets;
	std::vector<boost::filesystem::path> radon_matrix;
	boost::filesystem::path radon_matrix_second_factor;
	boost::filesystem::path rhs_file;
	int seed;
	boost::filesystem::path solution_file;
	boost::filesystem::path solution_image_file;

	// secondary values
};

#endif /* COMPUTERTOMOGRAPHYOPTIONS_HPP_ */
