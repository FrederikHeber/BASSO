// A simple program that computes the square root of a number
#include <iomanip>
#include <iostream>
#include <fstream>
#include "BassoConfig.h"

#include <boost/filesystem/path.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/program_options.hpp>

#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/LandweberMinimizer.hpp"
#include "Minimizations/SequentialSubspaceMinimizer.hpp"

namespace po = boost::program_options;
namespace logging = boost::log;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace keywords = boost::log::keywords;

BOOST_LOG_ATTRIBUTE_KEYWORD(timestamp, "TimeStamp", boost::posix_time::ptime);

int main (int argc, char *argv[])
{
	// set up command-line-parameters
	po::options_description desc("Allowed options");
	desc.add_options()
	        ("help", "produce help message")
	        ("verbose", po::value<unsigned int>(), "set the amount of verbosity")
	        ("C", po::value<double>(), "set the value for C")
	        ("normx", po::value<double>(), "set the norm of the space X")
	        ("normy", po::value<double>(), "set the norm of the space Y")
	        ("powery", po::value<double>(), "set the power type of the duality mapping's weight of the space Y")
	        ("delta", po::value<double>(), "set the amount of noise")
	        ("matrix", po::value< boost::filesystem::path >(),
	        		"set the forward operator matrix file")
	        ("rhs", po::value< boost::filesystem::path >(),
	        		"set the vector file of the right-hand side")
	        ("maxiter", po::value<unsigned int>(), "set the maximum amount of iterations")
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


	unsigned int verbose = 0;
	if (vm.count("verbose")) {
		verbose = vm["verbose"].as<unsigned int>();
		if (verbose > 0)
			std::cout << "Verbose set to " << verbose << std::endl;
	}
	switch (verbose) {
	default:
	case 0:
		logging::core::get()->set_filter
		(
				logging::trivial::severity >= logging::trivial::info
		);
		break;
	case 1:
		logging::core::get()->set_filter
		(
				logging::trivial::severity >= logging::trivial::debug
		);
		break;
	case 2:
		logging::core::get()->set_filter
		(
				logging::trivial::severity >= logging::trivial::trace
		);
		break;
	}
	// change log format
	logging::add_console_log(std::cout,
			keywords::format = (
				expr::stream
					<< expr::format_date_time(timestamp, "%Y-%m-%d %H:%M:%S")
					//<< expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
					<< ": <" << logging::trivial::severity
					<< "> " << expr::smessage
					)
			//keywords::format = "[%TimeStamp%]: [%Severity%] %Message%"
			);
	// register attributes
	logging::add_common_attributes();

	// parse options
	double C;
	if (vm.count("C")) {
		C = vm["C"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "C was set to " << C << ".\n";
	} else {
		C = 0.9;
	}
	double normx;
	if (vm.count("normx")) {
		normx = vm["normx"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Norm of X was set to " << normx << ".\n";
	}
	double normy;
	if (vm.count("normy")) {
		normy = vm["normy"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Norm of Y was set to " << normy << ".\n";
	}
	double powery;
	if (vm.count("powery")) {
		powery = vm["powery"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Power of duality maping in Y was set to " << powery << ".\n";
	}
	double delta;
	if (vm.count("delta")) {
		delta = vm["delta"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Magnitude of noise was set to " << delta << ".\n";
	}
	boost::filesystem::path matrix_file;
	if (vm.count("matrix")) {
		matrix_file = vm["matrix"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of matrix was set to " << matrix_file << ".\n";
	}
	boost::filesystem::path rhs_file;
	if (vm.count("rhs")) {
		rhs_file = vm["rhs"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of vector was set to " << rhs_file << ".\n";
	}
	unsigned int maxiter;
	if (vm.count("maxiter")) {
		maxiter = vm["maxiter"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Maximum iterations was set to " << maxiter << ".\n";
	} else {
		// set default value
		maxiter = 50;
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
	Eigen::VectorXd x0(rhs.innerSize());
	x0.setZero();
	std::cout << "Starting at x0 = " << x0.transpose() << std::endl;

	// call minimizer
//	SequentialSubspaceMinimizer minimizer(
//			normx,
//			normy,
//			powery,
//			delta,
//			maxiter);
//	SequentialSubspaceMinimizer::ReturnValues result =
//			minimizer(
//					x0,
//					matrix,
//					rhs);
	LandweberMinimizer minimizer(
			normx,
			normy,
			powery,
			delta,
			C,
			maxiter);
	LandweberMinimizer::ReturnValues result =
			minimizer(
					x0,
					matrix,
					rhs);

	// give result
	std::cout << "Solution is " << std::scientific << std::setprecision(8) << result.solution.transpose() << ","
		<< " found after " << result.NumberOuterIterations
		<< " with residual error of " << result.residuum << std::endl;

	// exit
	return 0;
}
