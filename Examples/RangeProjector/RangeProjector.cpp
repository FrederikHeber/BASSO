/*
 * RangeProjector.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <string>

#include <boost/assign.hpp>

#include "RangeProjector/Options/RangeProjectorOptions.hpp"
#include "RangeProjector/RangeProjector/RangeProjectionSolver.hpp"

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Solvers/SolverFactory/SolverFactory.hpp"

using namespace boost::assign;

int main (int argc, char *argv[])
{
	// starting timing
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	RangeProjectorOptions opts;
	opts.init();

	// parse options
	opts.parse(argc, argv);

	if (opts.showHelpConditions(argv[0]))
		return 1;

	// set verbosity level
	opts.setVerbosity();

	if (!opts.checkSensibility())
		return 255;
	opts.setSecondaryValues();

	// parse matrix and vector files into instances
	Eigen::MatrixXd matrix;
	Eigen::VectorXd rhs;
	{
		using namespace MatrixIO;

		{
			std::ifstream ist(opts.matrix_file.string().c_str());
			if (ist.good())
				try {
					ist >> matrix;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse matrix from " << opts.matrix_file.string() << std::endl;
					return 255;
				}
			else {
				std::cerr << "Failed to open " << opts.matrix_file.string() << std::endl;
				return 255;
			}

		}
		{
			std::ifstream ist(opts.rhs_file.string().c_str());
			if (ist.good())
				try {
					ist >> rhs;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse rhs from " << opts.rhs_file.string() << std::endl;
					return 255;
				}
			else {
				std::cerr << "Failed to open " << opts.rhs_file.string() << std::endl;
				return 255;
			}
		}
	}
	// print parsed matrix and vector if small or high verbosity requested
	if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
		BOOST_LOG_TRIVIAL(trace)
			<< "We solve for Ax = y with A = "
			<< matrix << " and y = "
			<< rhs.transpose() << std::endl;
	} else {
		BOOST_LOG_TRIVIAL(info)
			<< "We solve for Ax = y with A = "
			<< matrix << " and y = "
			<< rhs.transpose() << std::endl;
	}

	// create database
	Database_ptr_t database =
			SolverFactory::createDatabase(opts);

	// create projector instance
	RangeProjectionSolver projector(database,opts);

	// and project
	Eigen::VectorXd solution;
	projector(matrix, rhs, solution);

	// writing solution
	{
		using namespace MatrixIO;
		if (!opts.solution_file.string().empty()) {
			std::ofstream ost(opts.solution_file.string().c_str());
			if (ost.good())
				try {
					ost << std::setprecision(10) << solution;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write solution to file.\n";
				}
			else {
				std::cerr << "Failed to open " << opts.solution_file.string() << std::endl;
				return 255;
			}
		} else {
			std::cout << "No solution file given." << std::endl;
		}
	}

	boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	BOOST_LOG_TRIVIAL(info) << "The operation took "
			<< boost::chrono::duration<double>(timing_end - timing_start)
			<< ".";

	// exit
	return 0;
}



