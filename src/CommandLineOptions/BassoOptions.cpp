/*
 * BassoOptions.cpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include "BassoOptions.hpp"

#include <boost/filesystem.hpp>
#include <iostream>

#include "Log/Logging.hpp"

namespace po = boost::program_options;

BassoOptions::BassoOptions() :
			maxiter(50),
			maxwalltime(0.)
{}

void BassoOptions::internal_init()
{
	desc.add_options()
			("compare-against", po::value< boost::filesystem::path >(),
					"set the file name of the solution to compare against (BregmanDistance)")
	        ("matrix", po::value< boost::filesystem::path >(),
	        		"set the forward operator matrix file")
			("maxiter", po::value<unsigned int>(),
					"set the maximum amount of iterations")
			("max-walltime", po::value<double>(),
					"set the maximum time the algorithm may use")
	        ("rhs", po::value< boost::filesystem::path >(),
	        		"set the vector file of the right-hand side")
			("solution", po::value< boost::filesystem::path >(),
					"set the file name to write solution vector to, i.e. x in y = A*x")
			("solution-image", po::value< boost::filesystem::path >(),
					"set the file name to write image of solution vector to, i.e. A*x")
	        ;
}

void BassoOptions::internal_parse()
{
	if (vm.count("compare-against")) {
		comparison_file = vm["compare-against"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Parsing true solution vector from " << comparison_file;
	}

	if (vm.count("matrix")) {
		matrix_file = vm["matrix"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of matrix was set to " << matrix_file;
	}

	if (vm.count("maxiter")) {
		maxiter = vm["maxiter"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Maximum iterations was set to " << maxiter;
	}

	if (vm.count("max-walltime")) {
		maxwalltime = vm["max-walltime"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Maximum Walltime was set to " << maxwalltime;
	}

	if (vm.count("rhs")) {
		rhs_file = vm["rhs"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of vector was set to " << rhs_file;
	}

	if (vm.count("solution")) {
		solution_file = vm["solution"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Writing solution vector to " << solution_file;
	} else {
		BOOST_LOG_TRIVIAL(debug) << "No solution file given.";
	}

	if (vm.count("solution-image")) {
		solution_image_file = vm["solution-image"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Writing image of solution vector to " << solution_image_file;
	} else {
		BOOST_LOG_TRIVIAL(debug) << "No image of solution file given.";
	}
}

bool BassoOptions::internal_help_conditions() const
{
	if ((!vm.count("matrix")) || (!boost::filesystem::exists(matrix_file))) {
		BOOST_LOG_TRIVIAL(error)
				<< "Matrix file not set or not present.";
		return true;

	}

	if (!vm.count("normx")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Norm of space X normx not set";
		return true;
	}

	if (!vm.count("normy")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Norm of space Y normy not set";
		return true;
	}

	if ((!vm.count("rhs")) || (!boost::filesystem::exists(rhs_file))) {
		BOOST_LOG_TRIVIAL(error)
				<< "Right-hand side file not set or not present.";
		return true;

	}

	return false;
}

bool BassoOptions::internal_checkSensibility() const
{

	if ((vm.count("maxiter")) && (vm.count("max-walltime"))) {
		BOOST_LOG_TRIVIAL(error)
			<< "You have specified both max-iter and max-walltime.";
		return false;
	}

	return true;

}

void BassoOptions::internal_setSecondaryValues()
{}