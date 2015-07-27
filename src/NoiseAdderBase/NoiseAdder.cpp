/*
 * NoiseAdder.cpp
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

#include "NoiseAdderBase/Options/NoiseAdderOptions.hpp"

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "SolutionFactory/SolutionFactory.hpp"

using namespace boost::assign;

int main (int argc, char *argv[])
{
	// starting timing
	boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	NoiseAdderOptions opts;
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
	Eigen::MatrixXd input_element;
	{
		using namespace MatrixIO;

		{
			std::ifstream ist(opts.input_file.string().c_str());
			if (ist.good())
				try {
					ist >> input_element;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse rhs from " << opts.input_file.string() << std::endl;
					return 255;
				}
			else {
				std::cerr << "Failed to open " << opts.input_file.string() << std::endl;
				return 255;
			}
		}
	}

	// add random noise to matrix
	Eigen::MatrixXd output_element(input_element.rows(), input_element.cols());
	output_element.setRandom();
	if (opts.relativelevel) {
		const double maximum = input_element.maxCoeff();
		output_element *= maximum*opts.noiselevel;
	} else
		output_element *= opts.noiselevel;
	output_element += input_element;

	// writing solution
	{
		using namespace MatrixIO;
		if (!opts.output_file.string().empty()) {
			std::ofstream ost(opts.output_file.string().c_str());
			if (ost.good())
				try {
					ost << std::setprecision(10) << output_element;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write solution to file.\n";
				}
			else {
				std::cerr << "Failed to open " << opts.output_file.string() << std::endl;
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



