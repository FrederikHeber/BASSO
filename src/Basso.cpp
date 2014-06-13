
#include "BassoConfig.h"

// A simple program that computes the square root of a number
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/GeneralMinimizer.hpp"
#include "Minimizations/LandweberMinimizer.hpp"
#include "Minimizations/MinimizerFactory.hpp"
#include "Minimizations/SequentialSubspaceMinimizer.hpp"

namespace po = boost::program_options;

int main (int argc, char *argv[])
{
	// set up command-line-parameters
	po::options_description desc("Allowed options");
	desc.add_options()
			("algorithm", po::value<std::string>(), "set the used iteration algorithm")
	        ("C", po::value<double>(), "set the value for C (landweber)")
	        ("delta", po::value<double>(), "set the amount of noise")
	        ("help", "produce help message")
	        ("iteration-file", po::value< boost::filesystem::path >(),
	        		"set the filename to write information on iteration in sqlite format")
	        ("matrix", po::value< boost::filesystem::path >(),
	        		"set the forward operator matrix file")
			("maxiter", po::value<unsigned int>(), "set the maximum amount of iterations")
	        ("normx", po::value<double>(), "set the norm of the space X")
	        ("normy", po::value<double>(), "set the norm of the space Y")
	        ("number-directions", po::value<unsigned int>(), "set the number of search directions (SSO)")
	        ("output-steps", po::value<unsigned int>(), "output solution each ... steps")
	        ("powerx", po::value<double>(), "set the power type of the duality mapping's weight of the space X")
	        ("powery", po::value<double>(), "set the power type of the duality mapping's weight of the space Y")
	        ("rhs", po::value< boost::filesystem::path >(),
	        		"set the vector file of the right-hand side")
			("solution", po::value< boost::filesystem::path >(),
					"set the file name to write solution vector to")
			("compare-against", po::value< boost::filesystem::path >(),
					"set the file name of the solution to compare against (BregmanDistance)")
			("tau", po::value<double>(), "set the value for tau (SSO)")
	        ("verbose", po::value<unsigned int>(), "set the amount of verbosity")
	        ;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if ((vm.count("help"))
			|| (!vm.count("normx")) || (!vm.count("normy"))
			|| (!vm.count("delta"))
			|| (!vm.count("matrix")) || (!vm.count("rhs"))) {
		std::cout << argv[0] << " version "
				<< Basso_VERSION_MAJOR << "."
				<< Basso_VERSION_MINOR << std::endl;
	    std::cout << desc << "\n";
	    return 1;
	}

	// set verbosity level
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
	startLogging();

	// parse options
	MinimizerFactory::InstanceType type = MinimizerFactory::MAX_InstanceType;
	{
		// get desired algorithm
		std::string algorithm_name;
		if (vm.count("algorithm")) {
			algorithm_name = vm["algorithm"].as<std::string>();
		} else {
			algorithm_name =
					MinimizerFactory::TypeNames[MinimizerFactory::landweber];
		}
		// check whether algorithm_name states valid type
		if (!MinimizerFactory::isValidTypeName(algorithm_name)) {
			std::cerr << "Unknown algorithm specified by "
					<< algorithm_name << std::endl;
			return 255;
		}
		BOOST_LOG_TRIVIAL(debug)
			<< "algorithm was set to " << algorithm_name << "\n";
		type = MinimizerFactory::getTypeForName(algorithm_name);
	}
	double C;
	if (vm.count("C")) {
		C = vm["C"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "C was set to " << C << "\n";
	} else {
		C = 0.9;
	}
	double delta;
	boost::filesystem::path comparison_file;
	if (vm.count("compare-against")) {
		comparison_file = vm["compare-against"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Parsing true solution vector from " << comparison_file << ".\n";
	}
	if (vm.count("delta")) {
		delta = vm["delta"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Magnitude of noise was set to " << delta << "\n";
	}

	boost::filesystem::path iteration_file;
	if (vm.count("iteration-file")) {
		iteration_file = vm["iteration-file"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of iteration-file was set to " << iteration_file << "\n";
	}
	boost::filesystem::path matrix_file;
	if (vm.count("matrix")) {
		matrix_file = vm["matrix"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of matrix was set to " << matrix_file << "\n";
	}
	unsigned int maxiter;
	if (vm.count("maxiter")) {
		maxiter = vm["maxiter"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Maximum iterations was set to " << maxiter << "\n";
	} else {
		// set default value
		maxiter = 50;
	}
	double normx;
	if (vm.count("normx")) {
		normx = vm["normx"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Norm of X was set to " << normx << "\n";
	}
	double normy;
	if (vm.count("normy")) {
		normy = vm["normy"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Norm of Y was set to " << normy << "\n";
	}
	unsigned int N;
	if (vm.count("number-directions")) {
		N = vm["number-directions"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Number of search directions was set to " << N << "\n";
	} else {
		// set default value
		N = 2;
	}
	unsigned int outputsteps;
	if (vm.count("output-steps")) {
		outputsteps = vm["output-steps"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Output steps was set to " << outputsteps << "\n";
	} else {
		// set default value
		outputsteps = 0;
	}
	double powerx;
	if (vm.count("powerx")) {
		powerx = vm["powerx"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Power of duality maping in X was set to " << powerx << "\n";
	} else {
		BOOST_LOG_TRIVIAL(debug)
			<< "Using normx as powerx." << "\n.";
		powerx = normx;
	}
	double powery;
	if (vm.count("powery")) {
		powery = vm["powery"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Power of duality maping in Y was set to " << powery << "\n";
	} else {
		BOOST_LOG_TRIVIAL(debug)
			<< "Using normy as powery." << "\n.";
		powery = normy;
	}
	boost::filesystem::path rhs_file;
	if (vm.count("rhs")) {
		rhs_file = vm["rhs"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of vector was set to " << rhs_file << "\n";
	}
	double tau;
	if (vm.count("tau")) {
		tau = vm["tau"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "tau was set to " << tau << "\n";
	} else {
		tau = 1.1;
	}

	// parse matrix and vector files into instances
	Eigen::MatrixXd matrix;
	Eigen::VectorXd rhs;
	Eigen::VectorXd solution;
	{
		using namespace MatrixIO;

		{
			std::ifstream ist(matrix_file.string().c_str());
			if (ist.good())
				try {
					ist >> matrix;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse matrix from " << matrix_file.string() << std::endl;
					return 255;
				}
			else {
				std::cerr << "Failed to open " << matrix_file.string() << std::endl;
				return 255;
			}

		}
		{
			std::ifstream ist(rhs_file.string().c_str());
			if (ist.good())
				try {
					ist >> rhs;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse rhs from " << rhs_file.string() << std::endl;
					return 255;
				}
			else {
				std::cerr << "Failed to open " << rhs_file.string() << std::endl;
				return 255;
			}
		}
		{
			// just try to parse solution, we ignore if file does not exist
			std::ifstream ist(comparison_file.string().c_str());
			if (ist.good())
				try {
					ist >> solution;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully parse solution from " << comparison_file.string() << std::endl;
					return 255;
				}
			else {
				if (comparison_file.string().empty())
					BOOST_LOG_TRIVIAL(debug)
							<< "No solution file was given.";
				else {
					BOOST_LOG_TRIVIAL(error)
						<< "Could not parse solution from " << comparison_file.string();
					return 255;
				}
				solution = Eigen::VectorXd(matrix.outerSize());
				solution.setZero();
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

	// prepare start value
	Eigen::VectorXd x0(matrix.outerSize());
	x0.setZero();
	if (x0.innerSize() < 10)
		std::cout << "Starting at x0 = " << x0.transpose() << std::endl;

	// call minimizer
	MinimizerFactory factory;
	Database database;
	if (!iteration_file.string().empty())
		database.setDatabaseFile(iteration_file.string());
	MinimizerFactory::instance_ptr_t minimizer =
		factory.getInstance(
			type,
			normx,
			normy,
			powerx,
			powery,
			delta,
			maxiter,
			solution,
			database,
			outputsteps);
	try {
		switch(type) {
		case MinimizerFactory::landweber:
			static_cast<LandweberMinimizer*>(minimizer.get())->setC(C);
			break;
		case MinimizerFactory::sequentialsubspace:
			static_cast<SequentialSubspaceMinimizer*>(minimizer.get())->setTau(tau);
			static_cast<SequentialSubspaceMinimizer*>(minimizer.get())->setN(N);
			break;
		default:
			std::cerr << "Unknown InstanceType"
				<< MinimizerFactory::getNameForType(type) << "." << std::endl;
			return 255;
			break;
		}
	} catch (MinimizationIllegalValue_exception &e) {
		std::cerr << "Illegal value for "
				<< *boost::get_error_info<MinimizationIllegalValue_name>(e)
				<< std::endl;
		return 255;
	}
	GeneralMinimizer::ReturnValues result =
			(*minimizer)(
					x0,
					matrix,
					rhs);

	// give result
	if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
		std::cout << "Solution after "
				<< result.NumberOuterIterations
				<< " with residual error of " << result.residuum
				<< std::endl;
	} else {
		std::cout << "Solution after " << result.NumberOuterIterations
			<< " with residual error of " << result.residuum
			<< " is " << std::scientific << std::setprecision(8)
			<< result.solution.transpose()
			<< std::endl;
	}

	// writing solution
	{
		boost::filesystem::path solution_file;
		if (vm.count("solution")) {
			solution_file = vm["solution"].as<boost::filesystem::path>();
			BOOST_LOG_TRIVIAL(debug)
				<< "Writing solution vector to " << solution_file << ".\n";

		}

		{
			using namespace MatrixIO;
			std::ofstream ost(solution_file.string().c_str());
			if (ost.good())
				try {
					ost << result.solution;
				} catch (MatrixIOStreamEnded_exception &e) {
					std::cerr << "Failed to fully write solution to file.\n";
				}
			else {
				std::cerr << "Failed to open " << solution_file.string() << std::endl;
				return 255;
			}

		}
	}

	// exit
	return 0;
}
