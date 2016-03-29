/*
 * NonnegativeMatrixFactorsOptions.hpp
 *
 *  Created on: Nov 28, 2015
 *      Author: heber
 */

#ifndef NONNEGATIVEMATRIXFACTORSOPTIONS_HPP_
#define NONNEGATIVEMATRIXFACTORSOPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>

#include "Options/Options.hpp"

class NonnegativeMatrixFactorsOptions : public Options
{
public:
	NonnegativeMatrixFactorsOptions();

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
	boost::filesystem::path database_file;
	boost::filesystem::path destination_first_factor;
	boost::filesystem::path destination_second_factor;
	boost::filesystem::path matrix;
	unsigned int truncation_dimension;

	// secondary values
};

#endif /* NONNEGATIVEMATRIXFACTORSOPTIONS_HPP_ */
