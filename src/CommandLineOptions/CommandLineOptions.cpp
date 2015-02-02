/*
 * CommandLineOptions.cpp
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include "CommandLineOptions.hpp"

#include <iostream>
#include <sstream>

#include "Log/Logging.hpp"
#include "Minimizations/Minimizers/Searchspace/SearchspaceFactory.hpp"

namespace po = boost::program_options;

CommandLineOptions::CommandLineOptions() :
	desc("Allowed options"),
	algorithm_name(
			MinimizerFactory::TypeNames[MinimizerFactory::landweber]),
	C(.9),
	calculateAngles(false),
	database_replace(false),
	delta(1e-4),
	enforceRandomMapping(false),
	inexactLinesearch(false),
	maxinneriter(0),
	minlib("gsl"),
	normx(2.),
	normy(2.),
	N(2),
	outputsteps(0),
	powerx(2.),
	powery(2.),
	regularization_parameter(0),
	searchspace_type(
			SearchspaceFactory::getName(
					SearchspaceFactory::LastNDirections)),
	stepwidth_type(1),
	tau(1.1),
	updatetype(LastNSearchDirections::RoundRobin),
	dualitytype(defaulttype),
	type(MinimizerFactory::MAX_InstanceType)
{}

CommandLineOptions::~CommandLineOptions()
{}

void CommandLineOptions::init()
{
	desc.add_options()
			("algorithm", po::value<std::string>(),
					"set the used iteration algorithm")
			("C", po::value<double>(),
					"set the value for C (landweber)")
			("calculateAngles", po::value<bool>(),
					"set whether to calculate angles between search directions (SSO)")
			("database-replace", po::value<bool>(),
					"whether to replace tuples in database file (true) or just add (default=false)")
	        ("delta", po::value<double>(),
	        		"set the amount of noise")
			("enforceRandomMapping", po::value<bool>(),
					"whether to enforce update index algorithm in SSO to be a random mapping")
	        ("help", "produce help message")
	        ("inexact-linesearch", po::value<bool>(),
	        		"set the line search to be inexact or (quasi) exact")
	        ("iteration-file", po::value< boost::filesystem::path >(),
	        		"set the filename to write information on iteration in sqlite format")
			("minimization-library", po::value<std::string>(),
					"set which minimization library to use (gsl,nlopt)")
			("max-inner-iterations", po::value<unsigned int>(),
					"set the maximum amount of inner iterations")
	        ("normx", po::value<double>(),
	        		"set the norm of the space X, (1 <= p < inf, 0 means infinity norm)")
	        ("normy", po::value<double>(),
	        		"set the norm of the space Y, (1 <= r < inf, 0 means infinity norm)")
	        ("number-directions", po::value<unsigned int>(),
	        		"set the number of search directions (SSO)")
	        ("output-steps", po::value<unsigned int>(), "output solution each ... steps")
	        ("powerx", po::value<double>(),
	        		"set the power type of the duality mapping's weight of the space X")
	        ("powery", po::value<double>(),
	        		"set the power type of the duality mapping's weight of the space Y")
	        ("regularization-parameter", po::value<double>(),
	        		"set the regularization parameter for the L1 norm (if normx is 1) (adaptive if set to zero)")
			("searchspace", po::value< std::string >(),
					"set the type of search directions used (SSO")
			("stepwidth-algorithm", po::value<unsigned int>(),
					"set which step width algorithm to use (Landweber)")
			("tuple-parameters", po::value< std::vector<std::string> >()->multitoken(),
					"set additional parameters to add to tables in iteration database for distinguishing tuples")
			("tau", po::value<double>(),
					"set the value for tau (SSO)")
			("update-algorithm", po::value<unsigned int>(),
					"sets the algorithm which search direction is updated for multiple ones (SSO)")
	        ("verbose", po::value<unsigned int>(),
	        		"set the amount of verbosity")
			("wolfe-constants", po::value< std::vector<double> >()->multitoken(),
					"set the two wolfe conditions for positivity and stronger than linear (SSO, inexact line search)")
			;

	internal_init();
}

void CommandLineOptions::parse(int argc, char **argv)
{
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	// get desired algorithm
	if (vm.count("algorithm")) {
		algorithm_name = vm["algorithm"].as<std::string>();
		BOOST_LOG_TRIVIAL(debug)
			<< "algorithm was set to " << algorithm_name;
	}

	if (vm.count("C")) {
		C = vm["C"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "C was set to " << C;
	}

	if (vm.count("calculateAngles")) {
		calculateAngles = vm["calculateAngles"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "CalculateAngles was set to "
			<< (calculateAngles ? "true" : "false")
			<< ".";
	}

	if (vm.count("delta")) {
		delta = vm["delta"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Magnitude of noise was set to " << delta;
	}

	if (vm.count("database-replace")) {
		database_replace = vm["database-replace"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Database replace was set to "
			<< (database_replace ? "replace" : "just add") << ".";
	}

	if (vm.count("enforceRandomMapping")) {
		enforceRandomMapping = vm["enforceRandomMapping"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "We do " << (enforceRandomMapping ? "" : "not")
			<< " enforce the update algorithm to be a random mapping.";
	}

	if (vm.count("inexact-linesearch")) {
		inexactLinesearch = vm["inexact-linesearch"].as<bool>();
		BOOST_LOG_TRIVIAL(debug)
			<< "We do " << (inexactLinesearch ? "not do" : "")
			<< " an exact line search.";
	}

	if (vm.count("iteration-file")) {
		iteration_file = vm["iteration-file"].as<boost::filesystem::path>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Filename of iteration-file was set to " << iteration_file;
	}


	if (vm.count("max-inner-iterations")) {
		maxinneriter = vm["max-inner-iterations"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Maximum number of inner iterations was set to " << maxinneriter;
	}

	if (vm.count("minimization-library")) {
		minlib = vm["minimization-library"].as<std::string>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Using minimization library " << minlib;
	}

	if (vm.count("normx")) {
		normx = vm["normx"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Norm of X was set to " << normx;
		if (normx == 0.) // translate infinity
			normx = std::numeric_limits<double>::infinity();
	}

	if (vm.count("normy")) {
		normy = vm["normy"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Norm of Y was set to " << normy;
		if (normy == 0.) // translate infinity
			normy = std::numeric_limits<double>::infinity();
	}

	if (vm.count("number-directions")) {
		N = vm["number-directions"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Number of search directions was set to " << N;
	}

	if (vm.count("output-steps")) {
		outputsteps = vm["output-steps"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Output steps was set to " << outputsteps;
	}

	if (vm.count("powerx")) {
		powerx = vm["powerx"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Power of duality maping in X was set to " << powerx;
	} else {
		BOOST_LOG_TRIVIAL(debug)
			<< "Using normx as powerx." << "\n.";
		powerx = normx;
	}

	if (vm.count("powery")) {
		powery = vm["powery"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Power of duality maping in Y was set to " << powery;
	} else {
		BOOST_LOG_TRIVIAL(debug)
			<< "Using normy as powery." << "\n.";
		powery = normy;
	}

	if (vm.count("regularization-parameter")) {
		regularization_parameter = vm["regularization-parameter"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Regularization parameter of regularized L1 norm set to "
			<< regularization_parameter << ".";
	} else {
		BOOST_LOG_TRIVIAL(debug)
			<< "No regularization parameter set.";
	}

	if (vm.count("searchspace")) {
		searchspace_type = vm["searchspace"].as<std::string>();
		if (SearchspaceFactory::isValidName(searchspace_type)) {
			SearchspaceFactory::setCurrentType(
					SearchspaceFactory::getType(searchspace_type)
			);
			BOOST_LOG_TRIVIAL(debug)
				<< "Setting search space type to " << searchspace_type;
		}
	}

	if (vm.count("stepwidth-algorithm")) {
		if (type != MinimizerFactory::landweber) {
			BOOST_LOG_TRIVIAL(warning)
					<< "stepwidth-algorithm is set, but not landweber chosen, ignoring";
		}
		stepwidth_type = vm["stepwidth-algorithm"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug) << "stepwidth-algorithm was set to "
			<< stepwidth_type << "\n";
	}

	if (vm.count("tau")) {
		tau = vm["tau"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "tau was set to " << tau;
	}

	if (vm.count("tuple-parameters")) {
		tuple_parameters = vm["tuple-parameters"].as< std::vector<std::string> >();
		std::stringstream output;
		std::copy(tuple_parameters.begin(), tuple_parameters.end(),
				std::ostream_iterator<std::string>(output, ","));
		BOOST_LOG_TRIVIAL(debug)
			<< "tuple parameters was set to " << output.str();
	}

	if (vm.count("update-algorithm")) {
		updatetype = (enum LastNSearchDirections::UpdateAlgorithmType)
				vm["update-algorithm"].as<unsigned int>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Searchspace Index update algorithm set to " << updatetype;
	}

	if (vm.count("wolfe-constants")) {
		wolfe_constants = vm["wolfe-constants"].as< std::vector<double> >();
		BOOST_LOG_TRIVIAL(debug)
			<< "Setting wolfe positivity constant to " << wolfe_constants[0];
		BOOST_LOG_TRIVIAL(debug)
			<< "Setting wolfe stronger than linear constant to "
			<< wolfe_constants[1];
	}

	internal_parse();
}

bool CommandLineOptions::furtherHelpConditions() const
{
	if (!vm.count("delta")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Noise level delta not set";
		return true;
	}

	return false;
}

bool CommandLineOptions::showHelpConditions(const char * const program_name) const
{
	if (vm.count("help")) {
		std::cout << program_name << " version "
				<< Basso_VERSION_MAJOR << "."
				<< Basso_VERSION_MINOR << std::endl;
	    std::cout << desc << "\n";
	    return true;
	} else if (vm.count("version")) {
		std::cout << program_name << " version "
				<< Basso_VERSION_MAJOR << "."
				<< Basso_VERSION_MINOR << std::endl;
		return true;
	} else if (furtherHelpConditions() || internal_help_conditions()) {
		std::cout << "There was an error parsing options, use"
				<< "\n\t" << program_name << " --help\n"
				<< "to learn more." << std::endl;
	    return true;
	} else
		return false;
}

void CommandLineOptions::setVerbosity() const
{
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
}

bool CommandLineOptions::checkSensibility() const
{
	if (((normx == 1.) && !vm.count("regularization-parameter"))
		|| ((normx != 1.) && vm.count("regularization-parameter"))) {
		BOOST_LOG_TRIVIAL(error)
				<< "Either regularization parameter set but not l1 norm "
				<< "specified or the other way round.";
		return false;
	}
	if (vm.count("regularization-parameter")) {
		if (((vm.count("stepwidth-algorithm")
				&& (type != MinimizerFactory::landweber)))
				|| ((!vm.count("stepwidth-algorithm")
						&& (type == MinimizerFactory::landweber)))) {
			BOOST_LOG_TRIVIAL(error)
					<< "Either stepwidth algorithm chosen but not Landweber "
					<< "or the other way round";
			return false;
		}
	}

	if (vm.count("tau")) {
		const std::string &resesop_name =
				MinimizerFactory::getNameForType(
							MinimizerFactory::sequentialsubspace_noise);
		if (algorithm_name != resesop_name) {
			BOOST_LOG_TRIVIAL(warning)
					<< "Tau is set, but not " << resesop_name << ", ignoring";
			return false;
		}
	}

	// check whether algorithm_name states valid type
	if (!MinimizerFactory::isValidTypeName(algorithm_name)) {
		BOOST_LOG_TRIVIAL(error)
				<< "Unknown algorithm specified by " << algorithm_name;
		return false;
	}

	if ((minlib != "gsl") && (minlib != "nlopt")) {
		BOOST_LOG_TRIVIAL(error)
			<< "Minimization library must be either 'gsl' or 'nlopt'.";
		return false;
	}

	if (!SearchspaceFactory::isValidName(searchspace_type)) {
		BOOST_LOG_TRIVIAL(error)
					<< "Search space type " << searchspace_type
					<< " is unknown to factory.";
		return false;
	}

	if (vm.count("update-algorithm"))
		if (!((updatetype >= LastNSearchDirections::RoundRobin)
				&& (updatetype < LastNSearchDirections::MAX_UpdateAlgorithmType))) {
			BOOST_LOG_TRIVIAL(error)
					<< "Illegal update type set.";
			return false;
		}

	if (vm.count("wolfe-constants"))
		if ((wolfe_constants.size() != 2)
				|| (wolfe_constants[0] <= 0)
				|| (wolfe_constants[1] <= 0)) {
			BOOST_LOG_TRIVIAL(error)
					<< "Illegal Wolfe constants given, must be two and positive.";
			return false;
		}

	return internal_checkSensibility();
}


void CommandLineOptions::setSecondaryValues()
{
	// set dualitytype
	if (regularization_parameter > 0.)
		dualitytype = regularizedl1norm;
	else
		dualitytype = defaulttype;

	type = MinimizerFactory::getTypeForName(algorithm_name);

	// set good maximum of inner iterations
	if (maxinneriter == 0)
		maxinneriter = N*100;
}
