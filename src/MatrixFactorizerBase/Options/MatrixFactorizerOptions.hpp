/*
 * MatrixFactorizerOptions.hpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#ifndef MATRIXFACTORIZEROPTIONS_HPP_
#define MATRIXFACTORIZEROPTIONS_HPP_

#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>

#include "Options/CommandLineOptions.hpp"

#include <boost/serialization/access.hpp>
#include "filesystem_path_serialization.hpp"

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
	boost::filesystem::path data_file;
	unsigned int inner_iterations;
	unsigned int max_loops;
	double residual_threshold;
	boost::filesystem::path solution_factor_one_file;
	boost::filesystem::path solution_factor_two_file;
	boost::filesystem::path solution_product_file;
	unsigned int sparse_dim;

private:
	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive & ar, const unsigned int version) const
	{
		// MatrixFactorizerOptions
		ar & data_file;
		ar & inner_iterations;
		ar & max_loops;
		ar & residual_threshold;
		ar & solution_factor_one_file;
		ar & solution_factor_two_file;
		ar & solution_product_file;
		ar & sparse_dim;

		// CommandLineOptions
		ar & algorithm_name;
		ar & C;
		ar & calculateAngles;
		ar & config_filename;
		ar & database_replace;
		ar & delta;
		ar & enforceRandomMapping;
		ar & inexactLinesearch;
		ar & iteration_file;
		ar & maxinneriter;
		ar & minlib;
		ar & type_spacex;
		ar & type_spacey;
		ar & px;
		ar & py;
		ar & N;
		ar & orthogonalization_type;
		ar & outputsteps;
		ar & powerx;
		ar & powery;
		ar & regularization_parameter;
		ar & searchspace_type;
		ar & stepwidth_type;
		ar & tau;
		ar & tuple_parameters;
		ar & updatetype;
		ar & verbose;
		ar & wolfe_constants;
		ar & type;
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int version)
	{
		// MatrixFactorizerOptions
		ar & data_file;
		ar & inner_iterations;
		ar & max_loops;
		ar & residual_threshold;
		ar & solution_factor_one_file;
		ar & solution_factor_two_file;
		ar & solution_product_file;
		ar & sparse_dim;

		// CommandLineOptions
		ar & algorithm_name;
		ar & C;
		ar & calculateAngles;
		ar & config_filename;
		ar & database_replace;
		ar & delta;
		ar & enforceRandomMapping;
		ar & inexactLinesearch;
		ar & iteration_file;
		ar & maxinneriter;
		ar & minlib;
		ar & type_spacex;
		ar & type_spacey;
		ar & px;
		ar & py;
		ar & N;
		ar & orthogonalization_type;
		ar & outputsteps;
		ar & powerx;
		ar & powery;
		ar & regularization_parameter;
		ar & searchspace_type;
		ar & stepwidth_type;
		ar & tau;
		ar & tuple_parameters;
		ar & updatetype;
		ar & verbose;
		ar & wolfe_constants;
		ar & type;

		// set verbosity after loading
		setVerbosity();
	}

	BOOST_SERIALIZATION_SPLIT_MEMBER()
};



#endif /* MATRIXFACTORIZEROPTIONS_HPP_ */
