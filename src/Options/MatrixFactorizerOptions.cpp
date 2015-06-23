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
		inner_iterations(10),
		max_loops(100),
		residual_threshold(0.),
		sparse_dim(1)
{}

void MatrixFactorizerOptions::internal_init()
{
	boost::program_options::options_description desc_matrixfactorizer("MatrixFactorizer Options");

	desc_matrixfactorizer.add_options()
			("data", po::value<boost::filesystem::path>(),
					"set the file name of the data matrix")
			("inner-iterations", po::value<unsigned int>(),
					"set the maximum number of iterations spent on either matrix factor before switching")
			("max-loops", po::value<unsigned int>(),
					"set the maximum number of loops iterating over each factor")
			("residual-threshold", po::value<double>(),
					"set the threshold of the matrix residual when to stop the iteration")
			("solution-product", po::value< boost::filesystem::path >(),
					"set the file name to write the product of the two solution factors to")
			("solution-first-factor", po::value< boost::filesystem::path >(),
					"set the file name to write the first matrix factor")
			("solution-second-factor", po::value< boost::filesystem::path >(),
					"set the file name to write the second matrix factor")
			("sparse-dimension", po::value<unsigned int>(), "set the inner dimension of the matrix product")
			;

	desc_all.add(desc_matrixfactorizer);
}

void MatrixFactorizerOptions::internal_parse()
{
	if (vm.count("data")) {
		data_file = vm["data"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Data file name is " << data_file.string();
	}

	if (vm.count("inner-iterations")) {
		inner_iterations = vm["inner-iterations"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Performing " << inner_iterations << " inner iterations per factor.";
	}

	if (vm.count("max-loops")) {
		max_loops = vm["max-loops"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Performing " << max_loops << " loop iterations over the two factors.";
	}

	if (vm.count("residual-threshold")) {
		residual_threshold = vm["residual-threshold"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Stopping when matrix residual is less than " << residual_threshold << ".";
	}

	if (vm.count("solution-product")) {
		solution_product_file = vm["solution-product"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Solution product file name is " << solution_product_file.string();
	}

	if (vm.count("solution-first-factor")) {
		solution_factor_one_file = vm["solution-first-factor"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "First Solution factor file name is " << solution_factor_one_file.string();
	}

	if (vm.count("solution-second-factor")) {
		solution_factor_two_file = vm["solution-second-factor"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Second Solution factor file name is " << solution_factor_two_file.string();
	}

	if (vm.count("sparse-dimension")) {
		sparse_dim = vm["sparse-dimension"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Sparse dimension set to " << sparse_dim;
	}

}

bool MatrixFactorizerOptions::internal_help_conditions() const
{
	if ((!vm.count("data")) || (!boost::filesystem::exists(data_file))) {
		BOOST_LOG_TRIVIAL(error)
				<< "Data file not set or not present.";
		return true;

	}
	if (!vm.count("sparse-dimension")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Sparse dimensionality not set";
		return true;
	}

	return false;
}

bool MatrixFactorizerOptions::internal_checkSensibility() const
{
	// We have to check N+1 directions for linear independence. This should
	// be possible at least judging from the the dimensionality of the
	// space, hence this requirement. Otherwise, one of the offsets to
	// hyperplane becomes zero and the direction will just go for the
	// minimum norm (i.e. 0). This causes the iteration to jump between
	// N+1 different values, none of them the correct solution
	if ((N > 1) && (N > sparse_dim)) {
		BOOST_LOG_TRIVIAL(error)
				<< "Number of search directions must not be greater than the sparse dimension.";
		return false;
	}

	if (residual_threshold < 0.) {
		BOOST_LOG_TRIVIAL(error)
				<< "residual-threshold must be greater than non-negative";
		return false;
	}

	return true;
}

void MatrixFactorizerOptions::internal_setSecondaryValues()
{
	if (residual_threshold == 0.) {
		residual_threshold = delta;
		BOOST_LOG_TRIVIAL(info)
			<< "Residual threshold not set, defaulting to delta of "
			<< residual_threshold;
	}
}

void MatrixFactorizerOptions::internal_store(std::ostream &_output) const
{
	_output << "# [MatrixFactorizer]" << std::endl;
	writeValue<boost::filesystem::path>(_output, vm,  "data");
	writeValue<unsigned int>(_output, vm,  "inner-iterations");
	writeValue<unsigned int>(_output, vm,  "max-loops");
	writeValue<double>(_output, vm,  "residual-threshold");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-product");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-first-factor");
	writeValue<boost::filesystem::path>(_output, vm,  "solution-second-factor");
	writeValue<unsigned int>(_output, vm,  "sparse-dimension");
}
