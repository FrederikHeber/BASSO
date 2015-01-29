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
		max_loops(100)
{}

void MatrixFactorizerOptions::internal_init()
{
	desc.add_options()
			("data", po::value<boost::filesystem::path>(),
					"set the file name of the data matrix")
			("inner-iterations", po::value<unsigned int>(),
					"set the maximum number of iterations spent on either matrix factor before switching")
			("max-loops", po::value<unsigned int>(),
					"set the maximum number of loops iterating over each factor")
			("solution-product", po::value< boost::filesystem::path >(),
					"set the file name to write the product of the two solution factors to")
			("solution-first-factor", po::value< boost::filesystem::path >(),
					"set the file name to write the first matrix factor")
			("solution-second-factor", po::value< boost::filesystem::path >(),
					"set the file name to write the second matrix factor")
			("sparse-dimension", po::value<unsigned int>(), "set the inner dimension of the matrix product")
			;
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
	// We have to check N+1 directions for linear independency. This should
	// be possible at least judging from the the dimensionality of the
	// space, hence this requirement. Otherwise, one of the offsets to
	// hyperplane becoms zero and the direction will just go for the
	// minimum norm (i.e. 0). This causes the iteration to jump between
	// N+1 different values, none of them the correct solution
	if ((N > 1) && (N >= sparse_dim)) {
		BOOST_LOG_TRIVIAL(error)
				<< "Number of search directions must not be greater than the sparse dimension.";
		return false;
	}
	return true;
}
