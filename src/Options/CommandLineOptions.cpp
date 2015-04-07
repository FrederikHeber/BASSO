/*
 *
 *  Created on: Jan 27, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include "CommandLineOptions.hpp"

#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Log/Logging.hpp"
#include "Minimizations/Minimizers/Searchspace/SearchspaceFactory.hpp"

namespace po = boost::program_options;

CommandLineOptions::CommandLineOptions() :
	desc_all("Configuration Options"),
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
	boost::program_options::options_description desc_algorithm("Algorithm options");
	boost::program_options::options_description desc_banachspace("Banach space options");
	boost::program_options::options_description desc_general("General options");
	boost::program_options::options_description desc_landweber("Landweber options");
	boost::program_options::options_description desc_sesop("SESOP options");

	desc_algorithm.add_options()
			("algorithm", po::value<std::string>(),
					"set the used iteration algorithm")
	        ("delta", po::value<double>(),
	        		"set the amount of noise")
	        ("iteration-file", po::value< boost::filesystem::path >(),
	        		"set the filename to write information on iteration in sqlite format")
			("minimization-library", po::value<std::string>(),
					"set which minimization library to use (gsl,nlopt)")
			("max-inner-iterations", po::value<unsigned int>(),
					"set the maximum amount of inner iterations")
			;
	desc_banachspace.add_options()
			("normx", po::value<double>(),
					"set the norm of the space X, (1 <= p < inf, 0 means infinity norm)")
			("normy", po::value<double>(),
					"set the norm of the space Y, (1 <= r < inf, 0 means infinity norm)")
			("powerx", po::value<double>(),
					"set the power type of the duality mapping's weight of the space X")
			("powery", po::value<double>(),
					"set the power type of the duality mapping's weight of the space Y")
			("regularization-parameter", po::value<double>(),
					"set the regularization parameter for the L1 norm (if normx is 1) (adaptive if set to zero)")
			;

	desc_general.add_options()
			("database-replace", po::value<bool>(),
					"whether to replace tuples in database file (true) or just add (default=false)")
			("help", "produce help message")
	        ("output-steps", po::value<unsigned int>(), "output solution each ... steps")
			("tuple-parameters", po::value< std::vector<std::string> >()->multitoken(),
					"set additional parameters to add to tables in iteration database for distinguishing tuples")
			("verbose", po::value<unsigned int>(),
					"set the amount of verbosity")
			;

	desc_landweber.add_options()
			("C", po::value<double>(),
					"set the value for C")
			("stepwidth-algorithm", po::value<unsigned int>(),
					"set which step width algorithm to use")
			;

	desc_sesop.add_options()
			("number-directions", po::value<unsigned int>(),
					"set the number of search directions")
			("calculateAngles", po::value<bool>(),
					"set whether to calculate angles between search directions")
			("enforceRandomMapping", po::value<bool>(),
					"whether to enforce update index algorithm to be a random mapping")
			("inexact-linesearch", po::value<bool>(),
					"set the line search to be inexact or (quasi) exact")
			("searchspace", po::value< std::string >(),
					"set the type of search directions used")
			("tau", po::value<double>(),
					"set the value for discrepancy parameter tau")
			("update-algorithm", po::value<unsigned int>(),
					"sets the algorithm which search direction is updated for multiple ones")
			("wolfe-constants", po::value< std::vector<double> >()->multitoken(),
					"set the two wolfe conditions for positivity and stronger than linear (inexact line search)")
			;

	desc_config.add_options()
			("config", po::value<boost::filesystem::path>(),
					"filename of a configuration file containing default option parameters. Note that other command-line values override those in this file.")
			;

	desc_all
		.add(desc_config)
		.add(desc_general)
		.add(desc_banachspace)
		.add(desc_algorithm)
		.add(desc_landweber)
		.add(desc_sesop);

	internal_init();
}

void CommandLineOptions::parse(int argc, char **argv)
{
	/// parse the usual command line options
	po::store(po::parse_command_line(argc, argv, desc_all), vm);
	po::notify(vm);

	/// additionally parse config file
	if (vm.count("config")) {
		config_filename = vm["config"].as<boost::filesystem::path>();
		std::ifstream config_file(config_filename.string().c_str());
		po::store(po::parse_config_file(config_file, desc_all), vm);
		po::notify(vm);
	}

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
	} else if (vm.count("normx")) {
		BOOST_LOG_TRIVIAL(debug)
			<< "Using normx as powerx." << "\n.";
		powerx = normx;
	}

	if (vm.count("powery")) {
		powery = vm["powery"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Power of duality maping in Y was set to " << powery;
	} else if (vm.count("normy")) {
		BOOST_LOG_TRIVIAL(debug)
			<< "Using normy as powery." << "\n.";
		powery = normy;
	}

	if (vm.count("regularization-parameter")) {
		regularization_parameter = vm["regularization-parameter"].as<double>();
		BOOST_LOG_TRIVIAL(debug)
			<< "Regularization parameter of regularized L1 norm set to "
			<< regularization_parameter << ".";
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

bool CommandLineOptions::showHelpConditions(const char * const program_name) const
{
	if (vm.count("help")) {
		std::cout << program_name << " version "
				<< Basso_VERSION_MAJOR << "."
				<< Basso_VERSION_MINOR << std::endl;
	    std::cout << desc_all << "\n";
	    return true;
	} else if (vm.count("version")) {
		std::cout << program_name << " version "
				<< Basso_VERSION_MAJOR << "."
				<< Basso_VERSION_MINOR << std::endl;
		return true;
	} else
		return false;
}

void CommandLineOptions::showHelpinErrorCase() const
{
	std::cout << "There was an error parsing options, use '--help' to learn more." << std::endl;
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

bool CommandLineOptions::checkSensibility_delta() const
{
	if (!vm.count("delta")) {
		BOOST_LOG_TRIVIAL(error)
				<< "Noise level delta not set";
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_regularizationparameter() const
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
	return true;
}

bool CommandLineOptions::checkSensibility_tau() const
{
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
	return true;
}

bool CommandLineOptions::checkSensibility_tuple_parameters() const
{
	if (vm.count("tuple-parameters")) {
		if (tuple_parameters.size() % 2 == 0)
			return true;
		else
			return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_algorithm() const
{
	// check whether algorithm_name states valid type
	if (!MinimizerFactory::isValidTypeName(algorithm_name)) {
		BOOST_LOG_TRIVIAL(error)
				<< "Unknown algorithm specified by " << algorithm_name;
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_minlib() const
{
	if ((minlib != "gsl") && (minlib != "nlopt")) {
		BOOST_LOG_TRIVIAL(error)
			<< "Minimization library must be either 'gsl' or 'nlopt'.";
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_searchspace() const
{
	if (!SearchspaceFactory::isValidName(searchspace_type)) {
		BOOST_LOG_TRIVIAL(error)
					<< "Search space type " << searchspace_type
					<< " is unknown to factory.";
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_updatealgorithm() const
{
	if (vm.count("update-algorithm"))
		if (!((updatetype >= LastNSearchDirections::RoundRobin)
				&& (updatetype < LastNSearchDirections::MAX_UpdateAlgorithmType))) {
			BOOST_LOG_TRIVIAL(error)
					<< "Illegal update type set.";
			return false;
		}
	return true;
}

bool CommandLineOptions::checkSensibility_wolfeconstants() const
{
	if (vm.count("wolfe-constants"))
		if ((wolfe_constants.size() != 2)
				|| (wolfe_constants[0] <= 0)
				|| (wolfe_constants[1] <= 0)) {
			BOOST_LOG_TRIVIAL(error)
					<< "Illegal Wolfe constants given, must be two and positive.";
			return false;
		}
	return true;
}

bool CommandLineOptions::checkSensibility() const
{
	bool status = true;
	status &= checkSensibility_delta();
	status &= checkSensibility_regularizationparameter();
	status &= checkSensibility_tau();
	status &= checkSensibility_tuple_parameters();
	status &= checkSensibility_algorithm();
	status &= checkSensibility_minlib();
	status &= checkSensibility_searchspace();
	status &= checkSensibility_updatealgorithm();
	status &= checkSensibility_wolfeconstants();
	status &= internal_checkSensibility();

	if (!status)
		showHelpinErrorCase();

	return status;
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

void CommandLineOptions::store(std::ostream &_output)
{
	_output << "# [General]" << std::endl;
	writeValue<bool>(_output, vm,  "database-replace");
	writeValue<unsigned int>(_output, vm,  "output-steps");
	if (tuple_parameters.size() != 0) {
		for (std::vector<std::string>::const_iterator iter = tuple_parameters.begin();
				iter != tuple_parameters.end();)
			_output << "\ttuple-parameters = " << *(iter++) << std::endl;
	}
	writeValue<unsigned int>(_output, vm,  "verbose");

	_output << "# [Banach Space]" << std::endl;
	writeValue<double>(_output, vm,  "normx");
	writeValue<double>(_output, vm,  "normy");
	writeValue<double>(_output, vm,  "powerx");
	writeValue<double>(_output, vm,  "powery");
	writeValue<double>(_output, vm,  "regularization-parameter");

	_output << "# [Algorithm]" << std::endl;
	writeValue<std::string>(_output, vm,  "algorithm");
	writeValue<double>(_output, vm,  "delta");
	writeValue<boost::filesystem::path>(_output, vm,  "iteration-file");
	writeValue<std::string>(_output, vm,  "minimization-library");
	writeValue<unsigned int>(_output, vm,  "max-inner-iterations");

	_output << "# [Landweber]" << std::endl;
	writeValue<double>(_output, vm,  "C");
	writeValue<unsigned int>(_output, vm,  "stepwidth-algorithm");

	_output << "# [SESOP]" << std::endl;
	writeValue<unsigned int>(_output, vm,  "number-directions");
	writeValue<bool>(_output, vm,  "calculateAngles");
	writeValue<bool>(_output, vm,  "enforceRandomMapping");
	writeValue<bool>(_output, vm,  "inexact-linesearch");
	writeValue<std::string>(_output, vm,  "searchspace");
	writeValue<double>(_output, vm,  "tau");
	writeValue<unsigned int>(_output, vm,  "update-algorithm");
	if (vm.count("wolfe-constants") != 0) {
		const std::vector<double> values =
				vm["wolfe-constants"].as< std::vector<double> >();
		for (std::vector<double>::const_iterator iter = values.begin();
				iter != values.end();)
			_output << "\twolfe-constants = " << *(iter++) << std::endl;
	}

	internal_store(_output);
}
