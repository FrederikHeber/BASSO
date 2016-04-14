/*
 * MatrixFactorizerOptions.hpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZEROPTIONS_HPP_
#define MATRIXFACTORIZEROPTIONS_HPP_

#include "BassoConfig.h"

#include <vector>
#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>

#include "Options/filesystem_path_serialization.hpp"
#include "Options/CommandLineOptions.hpp"

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
	boost::filesystem::path data_file;
	unsigned int max_loops;
	std::vector<std::string> overall_keys;
	bool DoParseFactors;
	unsigned int fix_factor;
	double projection_delta;
	double residual_threshold;
	boost::filesystem::path solution_factor_one_file;
	boost::filesystem::path solution_factor_two_file;
	boost::filesystem::path solution_difference_file;
	boost::filesystem::path solution_product_file;
	bool sparse;
	unsigned int sparse_dim;

private:
	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int version) const
	{
		// CommandLineOptions
		ar & boost::serialization::base_object<const CommandLineOptions>(*this);

		// MatrixFactorizerOptions
		ar & data_file;
		ar & DoParseFactors;
		ar & fix_factor;
		ar & max_loops;
		ar & overall_keys;
		ar & projection_delta;
		ar & residual_threshold;
		ar & solution_difference_file;
		ar & solution_factor_one_file;
		ar & solution_factor_two_file;
		ar & solution_product_file;
		ar & sparse;
		ar & sparse_dim;
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int version)
	{
		// CommandLineOptions
		ar & boost::serialization::base_object<CommandLineOptions>(*this);

		// MatrixFactorizerOptions
		ar & data_file;
		ar & DoParseFactors;
		ar & fix_factor;
		ar & max_loops;
		ar & overall_keys;
		ar & projection_delta;
		ar & residual_threshold;
		ar & solution_difference_file;
		ar & solution_factor_one_file;
		ar & solution_factor_two_file;
		ar & solution_product_file;
		ar & sparse;
		ar & sparse_dim;

		// set secondary values
		setSecondaryValues();
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER()
};



#endif /* MATRIXFACTORIZEROPTIONS_HPP_ */
