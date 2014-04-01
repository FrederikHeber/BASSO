// A simple program that computes the square root of a number
#include <iostream>
#include <fstream>
#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>

#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/SequentialSubspaceMinimizer.hpp"

namespace po = boost::program_options;

int main (int argc, char *argv[])
{
	po::options_description desc("Allowed options");
	desc.add_options()
	        ("help", "produce help message")
	        ("normx", po::value<double>(), "set the norm of the space X")
	        ("normy", po::value<double>(), "set the norm of the space Y")
	        ("powery", po::value<double>(), "set the power type of the duality mapping's weight of the space Y")
	        ("delta", po::value<double>(), "set the amount of noise")
	        ("matrix", po::value< boost::filesystem::path >(),
	        		"set the forward operator matrix file")
	        ("rhs", po::value< boost::filesystem::path >(),
	        		"set the vector file of the right-hand side")
	        ;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if ((vm.count("help"))
			|| (!vm.count("normx")) || (!vm.count("normy"))
			|| (!vm.count("powery")) || (!vm.count("delta"))
			|| (!vm.count("matrix")) || (!vm.count("rhs"))) {
		std::cout << argv[0] << " version "
				<< Basso_VERSION_MAJOR << "."
				<< Basso_VERSION_MINOR << std::endl;
	    std::cout << desc << "\n";
	    return 1;
	}

	// parse options
	double normx;
	if (vm.count("normx")) {
		normx = vm["normx"].as<double>();
		std::cout << "Norm of X was set to " << normx << ".\n";
	}
	double normy;
	if (vm.count("normy")) {
		normy = vm["normy"].as<double>();
		std::cout << "Norm of Y was set to " << normy << ".\n";
	}
	double powery;
	if (vm.count("powery")) {
		powery = vm["powery"].as<double>();
		std::cout << "Power of duality maping in Y was set to " << powery << ".\n";
	}
	double delta;
	if (vm.count("delta")) {
		delta = vm["delta"].as<double>();
		std::cout << "Magnitude of noise was set to " << delta << ".\n";
	}
	boost::filesystem::path matrix_file;
	if (vm.count("matrix")) {
		matrix_file = vm["matrix"].as<boost::filesystem::path>();
		std::cout << "Filename of matrix was set to " << matrix_file << ".\n";
	}
	boost::filesystem::path rhs_file;
	if (vm.count("rhs")) {
		rhs_file = vm["rhs"].as<boost::filesystem::path>();
		std::cout << "Filename of vector was set to " << rhs_file << ".\n";
	}

	// parse matrix and vector files into instances
	Eigen::MatrixXd matrix;
	Eigen::VectorXd rhs;
	{
		using namespace MatrixIO;

		{
			std::ifstream ist(matrix_file.string().c_str());
			if (ist.good())
				ist >> matrix;
		}
		{
			std::ifstream ist(rhs_file.string().c_str());
			if (ist.good())
				ist >> rhs;
		}
	}
	std::cout << "We solve for Ax = y with A = "
			<< matrix << " and y = "
			<< rhs << "." << std::endl;

	// prepare start value
	Eigen::VectorXd x0(1, rhs.outerSize());
	x0.setZero();

	// call minimizer
	SequentialSubspaceMinimizer minimizer;
	SequentialSubspaceMinimizer::ReturnValues result =
			minimizer(
					x0,
					normx,
					matrix,
					rhs,
					normy,
					powery,
					delta);

	// give result
	std::cout << "Solution is " << result.solution << ","
			<< " found after " << result.NumberOuterIterations
			<< " with residual error of " << result.residuum << std::endl;

	// exit
	return 0;
}
