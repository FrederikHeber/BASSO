/*
 * MatrixFactorizerOptions.cpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include "MatrixFactorizerOptions.hpp"

#include <boost/filesystem.hpp>
#include <iostream>

#include "Log/Logging.hpp"

namespace po = boost::program_options;

MatrixFactorizerOptions::MatrixFactorizerOptions() :
		max_loops(100),
		DoParseFactors(false),
		fix_factor(0),
		projection_delta(1e-8),
		residual_threshold(0.),
		sparse_dim(1)
{}

void MatrixFactorizerOptions::internal_init()
{
	boost::program_options::options_description desc_matrixfactorizer("MatrixFactorizer Options");

	desc_matrixfactorizer.add_options()
			("data", po::value<boost::filesystem::path>(),
					"set the file name of the data matrix")
			("fix-factor", po::value< unsigned int >(),
					"set whether to fix one of the factors during minimization")
			("max-loops", po::value<unsigned int>(),
					"set the maximum number of loops iterating over each factor")
			("overall-keys", po::value< std::vector<std::string> >()->multitoken(),
					"set the columns to pass from inner problem's overall tables")
			("parse-factors", po::value< bool >(),
					"set whether to also parse initial factors from specified files")
			("projection-delta", po::value< double >(),
					"set the tolerance threshold for the projection of the variable column onto the range of the fixed matrix factor")
			("residual-threshold", po::value<double>(),
					"set the threshold of the matrix residual when to stop the iteration")
			("solution-difference", po::value< boost::filesystem::path >(),
					"set the file name to write the difference of the data matrix and the product of the two solution factors to")
			("solution-first-factor", po::value< boost::filesystem::path >(),
					"set the file name to write the first matrix factor (or also read in case of parse-factors)")
			("solution-product", po::value< boost::filesystem::path >(),
					"set the file name to write the product of the two solution factors to")
			("solution-second-factor", po::value< boost::filesystem::path >(),
					"set the file name to write the second matrix factor (or also read in case of parse-factors)")
			("sparse", po::value< bool >(),
					"set whether data matrix is parsed in sparse format")
			("sparse-dimension", po::value<unsigned int>(), "set the inner dimension of the matrix product")
			;

	desc_all.add(desc_matrixfactorizer);
}

void MatrixFactorizerOptions::internal_parse()
{
	if (vm.count("data")) {
		data_file = vm["data"].as<boost::filesystem::path>();
		LOG(debug, "Data file name is " << data_file.string());
	}
	if (vm.count("fix-factor")) {
		fix_factor = vm["fix-factor"].as<unsigned int>();
		if (fix_factor > 0) {
			LOG(debug, "We fix factor #" << fix_factor);
		} else {
			LOG(debug, "None of the factors are fixed.");
		}
	}

	if (vm.count("max-loops")) {
		max_loops = vm["max-loops"].as<unsigned int>();
		LOG(debug, "Performing " << max_loops << " loop iterations over the two factors.");
	}

	if (vm.count("overall-keys")) {
		overall_keys = vm["overall-keys"].as< std::vector<std::string> >();
		std::stringstream output;
		std::copy(overall_keys.begin(), overall_keys.end(),
				std::ostream_iterator<std::string>(output, ","));
		LOG(debug, "Overall keys was set to " << output.str());
	} else {
		overall_keys.clear();
		LOG(debug, "We don't accumulate any iteration keys from solver's overall table.");
	}
	if (vm.count("parse-factors")) {
		DoParseFactors = vm["parse-factors"].as<bool>();
		LOG(debug, "For the starting matrices, we " << (DoParseFactors ? "do": "do not") << " parse the factors from given files.");
	}
	if (vm.count("projection-delta")) {
		projection_delta = vm["projection-delta"].as<double>();
		LOG(debug, "Setting projection's delta to " << projection_delta);
	} else
		LOG(debug, "Setting projection's delta to default value of " << projection_delta);

	if (vm.count("residual-threshold")) {
		residual_threshold = vm["residual-threshold"].as<double>();
		LOG(debug, "Stopping when matrix residual is less than " << residual_threshold << ".");
	}

	if (vm.count("solution-difference")) {
		solution_difference_file = vm["solution-difference"].as<boost::filesystem::path>();
		LOG(debug, "Solution difference file name is " << solution_difference_file.string());
	}

	if (vm.count("solution-first-factor")) {
		solution_factor_one_file = vm["solution-first-factor"].as<boost::filesystem::path>();
		LOG(debug, "First Solution factor file name is " << solution_factor_one_file.string());
	}

	if (vm.count("solution-product")) {
		solution_product_file = vm["solution-product"].as<boost::filesystem::path>();
		LOG(debug, "Solution product file name is " << solution_product_file.string());
	}

	if (vm.count("solution-second-factor")) {
		solution_factor_two_file = vm["solution-second-factor"].as<boost::filesystem::path>();
		LOG(debug, "Second Solution factor file name is " << solution_factor_two_file.string());
	}

	if (vm.count("sparse-dimension")) {
		sparse_dim = vm["sparse-dimension"].as<unsigned int>();
		LOG(debug, "Sparse dimension set to " << sparse_dim);
	}
}

bool MatrixFactorizerOptions::internal_checkSensibility() const
{
	if ((!vm.count("data")) || (!boost::filesystem::exists(data_file))) {
		LOG(error, "Data file not set or not present.");
		return false;
	}

	if ((vm.count("fix-factor"))
			&& (fix_factor > 2)) {
		LOG(error, "Either first(1), second(2), or no (0) factors can be fixed.");
		return false;
	}

	if (!vm.count("sparse-dimension")) {
		LOG(error, "Sparse dimensionality not set");
		return false;
	}

	if (vm.count("parse-factors") && DoParseFactors) {
		// check files are given and exist
		if (!vm.count("solution-first-factor")
				|| !boost::filesystem::exists(solution_factor_one_file)) {
			LOG(error, "First factor file not specified or non-existent.");
			return false;
		}
		if (!vm.count("solution-second-factor")
				|| !boost::filesystem::exists(solution_factor_two_file)) {
			LOG(error, "Second factor file not specified or non-existent.");
			return false;
		}
	}

	// We have to check N+1 directions for linear independence. This should
	// be possible at least judging from the the dimensionality of the
	// space, hence this requirement. Otherwise, one of the offsets to
	// hyperplane becomes zero and the direction will just go for the
	// minimum norm (i.e. 0). This causes the iteration to jump between
	// N+1 different values, none of them the correct solution
	if ((N > 1) && (N > sparse_dim)) {
		LOG(error, "Number of search directions must not be greater than the sparse dimension.");
		return false;
	}

	if (residual_threshold < 0.) {
		LOG(error, "residual-threshold must be greater than non-negative");
		return false;
	}

	return true;
}

void MatrixFactorizerOptions::internal_setSecondaryValues()
{
	if (residual_threshold == 0.) {
		residual_threshold = delta;
		LOG(info, "Residual threshold not set, defaulting to delta of " << residual_threshold);
	}
}

void MatrixFactorizerOptions::internal_store(std::ostream &_output) const
{
	_output << "# [MatrixFactorizer]" << std::endl;
	writeValue<boost::filesystem::path>(_output, vm,  "data");
	writeValue<unsigned int>(_output, vm,  "fix-factor");
	writeValue<unsigned int>(_output, vm,  "max-loops");
	if (overall_keys.size() != 0) {
		for (std::vector<std::string>::const_iterator iter = overall_keys.begin();
				iter != overall_keys.end();)
			_output << "\toverall-keys = " << *(iter++) << std::endl;
	}
	writeValue<bool>(_output, vm,  "parse-factors");
	writeValue<double>(_output, vm, "projection-delta");
	writeValue<double>(_output, vm,  "residual-threshold");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-difference");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-first-factor");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-product");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-second-factor");
	writeValue<bool>(_output, vm,  "sparse");
	writeValue<unsigned int>(_output, vm,  "sparse-dimension");
}
