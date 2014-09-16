
#include "BassoConfig.h"

// A simple program that computes the square root of a number
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>

#include "Database/Database.hpp"
#include "Log/Logging.hpp"
#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/InverseProblems/InverseProblemFactory.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/MinimizerFactory.hpp"
#include "Minimizations/Minimizers/GeneralMinimizer.hpp"
#include "Minimizations/Minimizers/LandweberMinimizer.hpp"
#include "Minimizations/Minimizers/SequentialSubspaceMinimizer.hpp"
#include "Minimizations/Minimizers/SequentialSubspaceMinimizerNoise.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

namespace po = boost::program_options;


//!> enumeration of all known dualities (Don't forget to add string literal to TypeNames)
enum DualityContainerType {
	defaulttype=0,
	regularizedl1norm=1,
	MAX_DualityContainerType
};

int main (int argc, char *argv[])
{
	// set up command-line-parameters
	po::options_description desc("Allowed options");
	desc.add_options()
			("algorithm", po::value<std::string>(), "set the used iteration algorithm")
	        ("C", po::value<double>(), "set the value for C (landweber)")
			("compare-against", po::value< boost::filesystem::path >(),
					"set the file name of the solution to compare against (BregmanDistance)")
	        ("delta", po::value<double>(), "set the amount of noise")
	        ("enforceRandomMapping", po::value<bool>(), "whether to enforce update index algorithm in SSO to be a random mapping")
	        ("help", "produce help message")
	        ("inexact-linesearch", po::value<bool>(), "set the line search to be inexact or (quasi) exact")
	        ("iteration-file", po::value< boost::filesystem::path >(),
	        		"set the filename to write information on iteration in sqlite format")
	        ("matrix", po::value< boost::filesystem::path >(),
	        		"set the forward operator matrix file")
			("maxiter", po::value<unsigned int>(), "set the maximum amount of iterations")
	        ("normx", po::value<double>(), "set the norm of the space X, (1 <= p < inf, 0 means infinity norm)")
	        ("normy", po::value<double>(), "set the norm of the space Y, (1 <= r < inf, 0 means infinity norm)")
	        ("number-directions", po::value<unsigned int>(), "set the number of search directions (SSO)")
	        ("output-steps", po::value<unsigned int>(), "output solution each ... steps")
	        ("powerx", po::value<double>(), "set the power type of the duality mapping's weight of the space X")
	        ("powery", po::value<double>(), "set the power type of the duality mapping's weight of the space Y")
	        ("regularized-norm", po::value<double>(), "set the regularization parameter for the L1 norm (if normx is 1), activated if unequal zero")
	        ("rhs", po::value< boost::filesystem::path >(),
	        		"set the vector file of the right-hand side")
			("solution", po::value< boost::filesystem::path >(),
					"set the file name to write solution vector to")
			("tau", po::value<double>(), "set the value for tau (SSO)")
	        ("update-algorithm", po::value<unsigned int>(), "sets the algorithm which search direction is updated for multiple ones (SSO)")
			("useOptimalStepwidth", po::value<bool>(), "set whether to use optimal calculated step width or theoretical one (Landweber)")
	        ("verbose", po::value<unsigned int>(), "set the amount of verbosity")
			("wolfe-constants", po::value< std::vector<double> >()->multitoken(), "set the two wolfe conditions for positivity and stronger than linear (SSO, inexact line search)")
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
			<< "algorithm was set to " << algorithm_name;
		type = MinimizerFactory::getTypeForName(algorithm_name);
	}
	double C;
	if (vm.count("C")) {
		C = vm["C"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "C was set to " << C;
	} else {
		C = 0.9;
	}
	double delta;
	boost::filesystem::path comparison_file;
	if (vm.count("compare-against")) {
		comparison_file = vm["compare-against"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Parsing true solution vector from " << comparison_file;
	}
	if (vm.count("delta")) {
		delta = vm["delta"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Magnitude of noise was set to " << delta;
	} else {
		delta = 1e-4;
	}
	bool enforceRandomMapping = false;
	if (vm.count("enforceRandomMapping")) {
		enforceRandomMapping = vm["enforceRandomMapping"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "We do " << (enforceRandomMapping ? "" : "not")
			<< " enforce the update algorithm to be a random mapping.";
	}
	bool inexactLinesearch = false;
	if (vm.count("inexact-linesearch")) {
		inexactLinesearch = vm["inexact-linesearch"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "We do " << (inexactLinesearch ? "not do" : "")
			<< " an exact line search.";
	}
	boost::filesystem::path iteration_file;
	if (vm.count("iteration-file")) {
		iteration_file = vm["iteration-file"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of iteration-file was set to " << iteration_file;
	}
	boost::filesystem::path matrix_file;
	if (vm.count("matrix")) {
		matrix_file = vm["matrix"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of matrix was set to " << matrix_file;
	}
	unsigned int maxiter;
	if (vm.count("maxiter")) {
		maxiter = vm["maxiter"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Maximum iterations was set to " << maxiter;
	} else {
		// set default value
		maxiter = 50;
	}
	double normx;
	if (vm.count("normx")) {
		normx = vm["normx"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Norm of X was set to " << normx;
	} else {
		normx = 2;
	}
	if (normx == 0.) // translate infinity
		normx = std::numeric_limits<double>::infinity();
	double normy;
	if (vm.count("normy")) {
		normy = vm["normy"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Norm of Y was set to " << normy;
	} else {
		normy = 2.;
	}
	if (normy == 0.) // translate infinity
		normy = std::numeric_limits<double>::infinity();
	unsigned int N;
	if (vm.count("number-directions")) {
		N = vm["number-directions"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Number of search directions was set to " << N;
	} else {
		// set default value
		N = 2;
	}
	unsigned int outputsteps;
	if (vm.count("output-steps")) {
		outputsteps = vm["output-steps"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Output steps was set to " << outputsteps;
	} else {
		// set default value
		outputsteps = 0;
	}
	double powerx;
	if (vm.count("powerx")) {
		powerx = vm["powerx"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Power of duality maping in X was set to " << powerx;
	} else {
		BOOST_LOG_TRIVIAL(debug)
			<< "Using normx as powerx." << "\n.";
		powerx = normx;
	}
	double powery;
	if (vm.count("powery")) {
		powery = vm["powery"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Power of duality maping in Y was set to " << powery;
	} else {
		BOOST_LOG_TRIVIAL(debug)
			<< "Using normy as powery." << "\n.";
		powery = normy;
	}
	DualityContainerType dualitytype =
			defaulttype;
	double regularization_parameter = 0.;
	if ((vm.count("regularized-norm"))
			&& (vm["regularized-norm"].as<double>() > 0.)) {
		dualitytype = regularizedl1norm;
		regularization_parameter = vm["regularized-norm"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Duality Mappings of space X set to regularized L1 norm"
			<< " with parameter " << regularization_parameter << ".";
	} else {
		BOOST_LOG_TRIVIAL(debug)
			<< "Duality Mappings of space X set to default.";
	}
	boost::filesystem::path rhs_file;
	if (vm.count("rhs")) {
		rhs_file = vm["rhs"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of vector was set to " << rhs_file;
	}
	double tau;
	if (vm.count("tau")) {
		if (type != MinimizerFactory::sequentialsubspace_noise) {
			BOOST_LOG_TRIVIAL(warning)
					<< "Tau is set, but not SSO_noise chosen, ignoring";
		}
		tau = vm["tau"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "tau was set to " << tau;
	} else {
		tau = 1.1;
	}
	enum SequentialSubspaceMinimizer::UpdateAlgorithmType updatetype =
			SequentialSubspaceMinimizer::RoundRobin;
	if (vm.count("update-algorithm")) {
		const unsigned int temptype = vm["update-algorithm"].as<unsigned int>();
		if ((temptype >= SequentialSubspaceMinimizer::RoundRobin)
				&& (temptype < SequentialSubspaceMinimizer::MAX_UpdateAlgorithmType))
			updatetype =
					(enum SequentialSubspaceMinimizer::UpdateAlgorithmType)temptype;
		else {
			BOOST_LOG_TRIVIAL(error)
					<< "Illegal update type set.";
		}
		BOOST_LOG_TRIVIAL(debug)
			<< "Searchspace Index update algorithm set to " << updatetype;
	}
	bool useOptimalStepwidth;
	if (vm.count("useOptimalStepwidth")) {
		if (type != MinimizerFactory::landweber) {
			BOOST_LOG_TRIVIAL(warning)
					<< "useOptimalStepwidth is set, but not landweber chosen, ignoring";
		}
		useOptimalStepwidth = vm["useOptimalStepwidth"].as<bool>();
		BOOST_LOG_TRIVIAL(debug) << "useOptimalStepwidth was set to "
			<< useOptimalStepwidth;
	} else {
		useOptimalStepwidth = true;
	}
	std::vector<double> wolfe_constants;
	if (vm.count("wolfe-constants")) {
		wolfe_constants = vm["wolfe-constants"].as< std::vector<double> >();
		if ((wolfe_constants.size() != 2)
				|| (wolfe_constants[0] <= 0)
				|| (wolfe_constants[1] <= 0)) {
			std::cerr << "Illegal Wolfe constants given, must be two and positive."
					<< std::endl;
		}
		BOOST_LOG_TRIVIAL(debug)
			<< "Setting wolfe positivity constant to " << wolfe_constants[0];
		BOOST_LOG_TRIVIAL(debug)
			<< "Setting wolfe stronger than linear constant to "
			<< wolfe_constants[1];
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

	// prepare inverse problem
	InverseProblem_ptr_t inverseproblem;
	switch (dualitytype) {
	case regularizedl1norm:
		inverseproblem = InverseProblemFactory::createRegularizedL1Instance(
				regularization_parameter,
				powerx,
				normy,
				powery,
				matrix,
				rhs);
		break;
	case defaulttype:
	default:
		inverseproblem = InverseProblemFactory::createLpInstance(
				normx,
				powerx,
				normy,
				powery,
				matrix,
				rhs);
		break;
	}

	// prepare true solution
	SpaceElement_ptr_t truesolution =
			inverseproblem->x->getSpace()->createElement();
	*truesolution = solution;


	// prepare start value
	SpaceElement_ptr_t x0 =
			inverseproblem->x->getSpace()->createElement();
	x0->setZero();
	if (x0->getSpace()->getDimension() < 10)
		std::cout << "Starting at x0 = " << x0 << std::endl;

	// call minimizer
	MinimizerFactory factory;
	Database database;
	if (!iteration_file.string().empty())
		database.setDatabaseFile(iteration_file.string());
	MinimizerFactory::instance_ptr_t minimizer;
	// set regularization parametter in case of regularizedl1norm
	minimizer =
		factory.createInstance(
			type,
			inverseproblem,
			delta,
			maxiter,
			database,
			outputsteps);

	// calculate initial dual solution
	SpaceElement_ptr_t dualx0 =
			(dualitytype == defaulttype) ?
			(*inverseproblem->x->getSpace()->getDualityMapping())(x0) :
			inverseproblem->x->getSpace()->getDualSpace()->createElement();

	try {
		switch(type) {
		case MinimizerFactory::landweber:
			static_cast<LandweberMinimizer*>(
					minimizer.get())->setC(C);
			static_cast<LandweberMinimizer*>(
					minimizer.get())->setuseOptimalStepwidth(
					useOptimalStepwidth);
			break;
		case MinimizerFactory::sequentialsubspace:
			static_cast<SequentialSubspaceMinimizer*>(minimizer.get())->setN(N);
			static_cast<SequentialSubspaceMinimizer*>(
					minimizer.get())->setupdateIndexAlgorithm(
							inverseproblem,
							updatetype);
			static_cast<SequentialSubspaceMinimizer*>(
					minimizer.get())->setEnforceRandomMapping(
							enforceRandomMapping);
			static_cast<SequentialSubspaceMinimizer*>(minimizer.get())->setInexactLinesearch(
					inexactLinesearch);
			if (wolfe_constants.size() == 2) {
				static_cast<SequentialSubspaceMinimizer*>(minimizer.get())->setWolfeConstants(
						wolfe_constants);
				// warning in case sanity check fails
				if (!inexactLinesearch)
					BOOST_LOG_TRIVIAL(warning)
						<< "Wolfe constants set although we do perform an exact line search.";
			}
			break;
		case MinimizerFactory::sequentialsubspace_noise:
			static_cast<SequentialSubspaceMinimizerNoise*>(
					minimizer.get())->setTau(tau);
			static_cast<SequentialSubspaceMinimizerNoise*>(
					minimizer.get())->setN(N);
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
					inverseproblem,
					x0,
					dualx0,
					truesolution);
	minimizer->resetState();

	// give result
	{
		Norm_ptr_t NormY = inverseproblem->y->getSpace()->getNorm();
		if ((matrix.innerSize() > 10) || (matrix.outerSize() > 10)) {
			std::cout << "Solution after "
					<< result.NumberOuterIterations
					<< " with relative residual of " << result.residuum/(*NormY)(rhs)
					<< std::endl;
		} else {
			std::cout << "Solution after " << result.NumberOuterIterations
				<< " with relative residual of " << result.residuum/(*NormY)(rhs)
				<< " is " << std::scientific << std::setprecision(8)
				<< inverseproblem->x
				<< std::endl;
		}
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
					ost << result.m_solution->getVectorRepresentation();
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
