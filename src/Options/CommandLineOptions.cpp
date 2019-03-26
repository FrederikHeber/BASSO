/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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

#include <boost/filesystem.hpp>

#include "Log/Logging.hpp"
#include "Minimizations/Functions/Minimizers/FunctionalMinimizerFactory.hpp"
#include "Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp"
#include "Minimizations/Minimizers/Searchspace/SearchspaceFactory.hpp"
#include "Minimizations/Minimizers/StepWidths/DetermineStepWidthFactory.hpp"
#include "Minimizations/Norms/NormFactory.hpp"

namespace po = boost::program_options;

CommandLineOptions::CommandLineOptions() :
	algorithm_name(
			MinimizerFactory::TypeNames[MinimizerFactory::landweber]),
	auxiliary_constraints(""),
	C(.9),
	calculateAngles(false),
	database_replace(false),
	delta(1e-4),
	enforceRandomMapping(false),
	everynthtuple(1),
	inexactLinesearch(false),
	maxinneriter(0),
	maxiter(50),
	max_sfp_loops(1),
	maxwalltime(0.),
	minlib("gsl"),
	type_spacex("lp"),
	type_spacey("lp"),
	px(2.),
	py(2.),
	N(1),
	orthogonalization_type(LastNSearchDirections::NoOrthogonalization),
	outputsteps(0),
	powerx(2.),
	powery(2.),
	regularization_parameter(0),
	searchspace_type(
			SearchspaceFactory::getName(
					SearchspaceFactory::LastNDirections)),
	stepwidth_type(DetermineStepWidthFactory::MinimizingResidual),
	tau(1.),
	tolerance_linesearch(1e-12),
	tolerance_spacex(1e-6),
	updatetype(LastNSearchDirections::RoundRobin),
	type(MinimizerFactory::MAX_InstanceType)
{}

CommandLineOptions::~CommandLineOptions()
{}

void CommandLineOptions::init()
{
	Options::init();

	boost::program_options::options_description desc_algorithm("Algorithm options");
	boost::program_options::options_description desc_banachspace("Banach space options");
	boost::program_options::options_description desc_general("General options");
	boost::program_options::options_description desc_landweber("Landweber options");
	boost::program_options::options_description desc_sesop("SESOP options");
	boost::program_options::options_description desc_sfp("Split feasibility options");

	desc_algorithm.add_options()
			("algorithm", po::value<std::string>(),
					"set the used iteration algorithm")
			("delta", po::value<double>(),
					"set (optionally) the tolerance value for performing dual mapping of elements in space Y.")
			("iteration-file", po::value< boost::filesystem::path >(),
	        		"set the filename to write information on iteration in sqlite format")
			("every-nth-tuple", po::value<unsigned int>(),
					"do not store every iteration information tuple in iteration file but only every nth, 0 deactivates per_iteration tuples.")
			("minimization-library", po::value<std::string>(),
					"set which minimization library to use (gsl,nlopt)")
			("max-inner-iterations", po::value<unsigned int>(),
					"set the maximum amount of inner iterations")
			("maxiter", po::value<unsigned int>(),
					"set the maximum amount of iterations")
			("max-walltime", po::value<double>(),
					"set the maximum time the algorithm may use")
			("stopping-criteria", po::value< std::string >(),
					"set (optionally) the desired stopping criteria, boolean operators allowed. This overrules any set stopping parameters such as maxiter, max-walltime, and delta.")
			("tolerance-space-x", po::value< double >(),
					"set (optionally) the tolerance value for performing dual mapping of elements in space X.")
			("tolerance-linesearch", po::value< double >(),
					"set (optionally) the tolerance value for line searches.")
			;
	desc_banachspace.add_options()
			("type-space-x", po::value<std::string>(),
					"sets the type of the norm for the source space X")
			("type-space-y", po::value<std::string>(),
					"sets the type of the norm for the destination space Y")
			("px", po::value<double>(),
					"set the norm of the space X, (1 <= p < inf, 0 means infinity norm)")
			("py", po::value<double>(),
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
			("output-steps", po::value<unsigned int>(),
					"output solution each ... steps")
			("tuple-parameters", po::value< std::vector<std::string> >()->multitoken(),
					"set additional parameters to add to tables in iteration database for distinguishing tuples")
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
			("orthogonal-directions", po::value<unsigned int>(),
					"set the way of orthogonalizing search directions (none(0), metric(1), bregman(2))")
			("searchspace", po::value< std::string >(),
					"set the type of search directions used")
			("tau", po::value<double>(),
					"set the value for discrepancy parameter tau, default is 1.1")
			("update-algorithm", po::value<unsigned int>(),
					"sets the algorithm which search direction is updated for multiple ones")
			("wolfe-constants", po::value< std::vector<double> >()->multitoken(),
					"set the two wolfe conditions for positivity and stronger than linear (inexact line search)")
			;
	desc_sfp.add_options()
			("auxiliary-constraints", po::value<std::string>(),
					"set any auxiliary constraints")
			("max-sfp-loops", po::value<unsigned int>(),
					"set the maximum number of loops for Split Feasibility Problems")
			;

	desc_all
		.add(desc_general)
		.add(desc_banachspace)
		.add(desc_algorithm)
		.add(desc_landweber)
		.add(desc_sesop)
		.add(desc_sfp);

	internal_init();
}

void CommandLineOptions::parse(int argc, char **argv)
{
	Options::parse(argc,argv);

	// get desired algorithm
	if (vm.count("algorithm")) {
		algorithm_name = vm["algorithm"].as<std::string>();
		LOG(debug, "algorithm was set to " << algorithm_name);
	}

	if (vm.count("auxiliary-constraints")) {
		auxiliary_constraints = vm["auxiliary-constraints"].as<std::string>();
		LOG(debug, "Auxiliary constraints are " << (auxiliary_constraints.empty() ? "None" : auxiliary_constraints));
	} else {
		LOG(debug, "Auxiliary constraints are None");
	}

	if (vm.count("C")) {
		C = vm["C"].as<double>();
		LOG(debug, "C was set to " << C);
	}

	if (vm.count("calculateAngles")) {
		calculateAngles = vm["calculateAngles"].as<bool>();
		LOG(debug, "CalculateAngles was set to "
			<< (calculateAngles ? "true" : "false") << ".");
	}

	if (vm.count("delta")) {
		delta = vm["delta"].as<double>();
		LOG(debug, "Magnitude of noise was set to " << delta);
	} else {
		LOG(debug, "delta set to default value of " << delta);
	}

	if (vm.count("database-replace")) {
		database_replace = vm["database-replace"].as<bool>();
		LOG(debug, "Database replace was set to " << (database_replace ? "replace" : "just add") << ".");
	}

	if (vm.count("enforceRandomMapping")) {
		enforceRandomMapping = vm["enforceRandomMapping"].as<bool>();
		LOG(debug, "We do " << (enforceRandomMapping ? "" : "not") << " enforce the update algorithm to be a random mapping.");
	}

	if (vm.count("every-nth-tuple")) {
		everynthtuple = vm["every-nth-tuple"].as<unsigned int>();
		if (everynthtuple != 0) {
			LOG(debug, "Every nth tuple was set to " << everynthtuple);
		} else {
			LOG(info, "Per iteration tuples are not added in database, only overall information.");
		}
	} else {
		LOG(debug, "Storing every tuple in database.");
	}

	if (vm.count("inexact-linesearch")) {
		inexactLinesearch = vm["inexact-linesearch"].as<bool>();
		LOG(debug, "We do " << (inexactLinesearch ? "not do" : "") << " an exact line search.");
	}

	if (vm.count("iteration-file")) {
		iteration_file = vm["iteration-file"].as<boost::filesystem::path>();
		LOG(debug, "Filename of iteration-file was set to " << iteration_file);
	}

	if (vm.count("maxiter")) {
		maxiter = vm["maxiter"].as<unsigned int>();
		LOG(debug, "Maximum iterations was set to " << maxiter);
	} else {
		LOG(debug, "Maximum iterations set to default value of " << maxiter);
	}

	if (vm.count("max-inner-iterations")) {
		maxinneriter = vm["max-inner-iterations"].as<unsigned int>();
		LOG(debug, "Maximum number of inner iterations was set to " << maxinneriter);
	}

	if (vm.count("max-sfp-loops")) {
		max_sfp_loops = vm["max-sfp-loops"].as<unsigned int>();
		LOG(debug, "Maximum number of SFP loops was set to " << max_sfp_loops);
	} else {
		LOG(debug, "Maximum of SFP loops set to default value of " << max_sfp_loops);
	}

	if (vm.count("max-walltime")) {
		maxwalltime = vm["max-walltime"].as<double>();
		LOG(debug, "Maximum Walltime was set to " << maxwalltime << " seconds.");
	} else {
		LOG(debug, "Maximum Walltime set to default value of " << boost::chrono::duration<double>(maxwalltime) << " seconds.");
	}

	if (vm.count("minimization-library")) {
		minlib = vm["minimization-library"].as<std::string>();
		LOG(debug, "Using minimization library " << minlib);
	}

	if (vm.count("type-space-x")) {
		type_spacex = vm["type-space-x"].as<std::string>();
		LOG(debug, "Type of space X was set to " << type_spacex);
	}

	if (vm.count("type-space-y")) {
		type_spacey = vm["type-space-y"].as<std::string>();
		LOG(debug, "Type of space Y was set to " << type_spacey);
	}

	if (vm.count("px")) {
		px = vm["px"].as<double>();
		LOG(debug, "P value of lp norm of X was set to " << px);
		if (px == 0.) // translate infinity
			px = std::numeric_limits<double>::infinity();
	}

	if (vm.count("py")) {
		py = vm["py"].as<double>();
		LOG(debug, "P value of lp norm of Y was set to " << py);
		if (py == 0.) // translate infinity
			py = std::numeric_limits<double>::infinity();
	}

	if (vm.count("number-directions")) {
		N = vm["number-directions"].as<unsigned int>();
		LOG(debug, "Number of search directions was set to " << N);
	}

	if (vm.count("orthogonal-directions")) {
		orthogonalization_type =
				(enum LastNSearchDirections::OrthogonalizationType)
						vm["orthogonal-directions"].as<unsigned int>();
		LOG(debug, "Orthogonalizing search directions: " << LastNSearchDirections::getOrthogonalizationTypeName(orthogonalization_type));
	}

	if (vm.count("output-steps")) {
		outputsteps = vm["output-steps"].as<unsigned int>();
		LOG(debug, "Output steps was set to " << outputsteps);
	}

	if (vm.count("powerx")) {
		powerx = vm["powerx"].as<double>();
		LOG(debug, "Power of duality maping in X was set to " << powerx);
	} else if (vm.count("px")) {
		powerx = px > 2. ? px : 2.;
		LOG(debug, "Using px as powerx if above 2: " << powerx << "\n.");
	}

	if (vm.count("powery")) {
		powery = vm["powery"].as<double>();
		LOG(debug, "Power of duality maping in Y was set to " << powery);
	} else if (vm.count("py")) {
		powery = py > 2. ? py : 2.;
		LOG(debug, "Using py as powery if above 2: " << powery << "\n.");
	}

	if (vm.count("regularization-parameter")) {
		regularization_parameter = vm["regularization-parameter"].as<double>();
		LOG(debug, "Regularization parameter of regularized L1 norm set to " << regularization_parameter << ".");
	}

	if (vm.count("searchspace")) {
		searchspace_type = vm["searchspace"].as<std::string>();
		LOG(debug, "Setting search space type to " << searchspace_type);
	}

	if (vm.count("stepwidth-algorithm")) {
		stepwidth_type = vm["stepwidth-algorithm"].as<unsigned int>();
		LOG(debug, "stepwidth-algorithm was set to " << stepwidth_type << "\n");
	}

	if (vm.count("stopping-criteria")) {
		stopping_criteria = vm["stopping-criteria"].as< std::string >();
		LOG(debug, "Stopping criteria was set to " << stopping_criteria);
	}

	if (vm.count("tau")) {
		tau = vm["tau"].as<double>();
		LOG(debug, "tau was set to " << tau);
	} else {
		LOG(debug, "tau set to default value of " << tau);
	}

	if (vm.count("tolerance-linesearch")) {
		tolerance_linesearch = vm["tolerance-linesearch"].as<double>();
		LOG(debug, "Tolerance for linesearches was set to " << tolerance_linesearch);
	} else {
		LOG(debug, "Tolerance for linesearches set to default value of " << tolerance_linesearch);
	}

	if (vm.count("tolerance-space-x")) {
		tolerance_spacex = vm["tolerance-space-x"].as<double>();
		LOG(debug, "Tolerance for space x mappings was set to " << tolerance_spacex);
	} else {
		LOG(debug, "Tolerance for space x mappings set to default value of " << tolerance_spacex);
	}

	if (vm.count("tuple-parameters")) {
		tuple_parameters = vm["tuple-parameters"].as< std::vector<std::string> >();
		std::stringstream output;
		std::copy(tuple_parameters.begin(), tuple_parameters.end(),
				std::ostream_iterator<std::string>(output, ","));
		LOG(debug, "tuple parameters was set to " << output.str());
	}

	if (vm.count("update-algorithm")) {
		updatetype = (enum LastNSearchDirections::UpdateAlgorithmType)
				vm["update-algorithm"].as<unsigned int>();
		LOG(debug, "Searchspace Index update algorithm set to " << updatetype);
	}

	if (vm.count("wolfe-constants")) {
		wolfe_constants = vm["wolfe-constants"].as< std::vector<double> >();
		LOG(debug, "Setting wolfe positivity constant to " << wolfe_constants[0]);
		LOG(debug, "Setting wolfe stronger than linear constant to " << wolfe_constants[1]);
	}

	internal_parse();
}

bool CommandLineOptions::checkSensibility_config() const
{
	if (vm.count("config")) {
		if (!boost::filesystem::exists(config_filename)) {
				LOG(error, "Specified config file " << config_filename << " does not exist.");
				return false;
		}
	}
	return true;
}

bool CommandLineOptions::checkSensibility_delta() const
{
	if (!vm.count("delta")) {
		LOG(error, "Noise level delta not set");
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_everynthtuple() const
{
	if (!vm.count("every-nth-tuple")) {
		if (vm.count("maxiter") && (maxiter > 5000)) {
			LOG(error, "Maxiter is very large. We recommend to not store all tuples, "
					<< "use every-nth-tuple. Or set explicitly to 1 if desired");
			return false;
		}
	}
	return true;
}

bool CommandLineOptions::checkSensibility_OrthogonalDirections() const
{
	if (vm.count("orthogonal-directions")) {
		if (!((orthogonalization_type >= LastNSearchDirections::NoOrthogonalization)
			&& (orthogonalization_type < LastNSearchDirections::MAX_Orthogonalization))) {
			LOG(error, "Illegal orthogonalization type specified.");
			return false;
		}
	}
	return true;
}

bool CommandLineOptions::checkSensibility_regularizationparameter() const
{
	if (type_spacex == "regularized_l1") {
		if (!vm.count("regularization-parameter")) {
			LOG(error, "No regularization parameter set.");
			return false;
		} else {
			if (vm.count("px") && (px != 1.)) {
				LOG(error, "Regularization parameter given and px set but not to l1 norm.");
				return false;
			}
		}
	}

	if (type_spacey == "regularized_l1") {
		if (!vm.count("regularization-parameter")) {
			LOG(error, "No regularization parameter set.");
			return false;
		} else {
			if (vm.count("py") && (py != 1.)) {
				LOG(error, "Regularization parameter given and py set but not to l1 norm.");
				return false;
			}
		}
	}

	return true;
}

bool CommandLineOptions::checkSensibility_tau() const
{
	if ((vm.count("tau")) && (tau < 1.)) {
		LOG(warning, "Tau is set, but value of " << tau << " is not greater equal to 1.");
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_tuple_parameters() const
{
	if (vm.count("tuple-parameters")) {
		if (tuple_parameters.size() % 2 == 0)
			return true;
		else {
			LOG(error, "Tuple parameters need to come in pairs.");
			return false;
		}
	}
	return true;
}

bool CommandLineOptions::checkSensibility_algorithm() const
{
	// check whether algorithm_name states valid type
	if (!MinimizerFactory::isValidTypeName(algorithm_name)) {
		LOG(error, "Unknown algorithm specified by " << algorithm_name);
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_minlib() const
{
	if ((minlib != "gsl") && (minlib != "nlopt")) {
		LOG(error, "Minimization library must be either 'gsl' or 'nlopt'.");
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_norms() const
{
	if ((!NormFactory::isValidType(type_spacex))
			|| (!NormFactory::isValidType(type_spacey))) {
		LOG(error, "Type of space of X or Y has invalid type.");
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_pvalues() const
{
	if ((type_spacex.substr(0,2) != "lp") &&
			((vm.count("px") != 0) && (px != 1.))) {
		LOG(error, "No lp type desired but px value unequal 1 given");
		return false;
	}
	if ((type_spacey.substr(0,2) != "lp") &&
			((vm.count("py") != 0) && (py != 1.))) {
		LOG(error, "No lp type desired but py value unequal 1 given");
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_searchspace() const
{
	if (!SearchspaceFactory::isValidName(searchspace_type)) {
		LOG(error, "Search space type " << searchspace_type << " is unknown to factory.");
		return false;
	}

	if (!FunctionalMinimizerFactory::isValidName(minlib)) {
		LOG(error, "Minimization library " << minlib << " is unknown to factory.");
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_stepwidth_algorithm() const
{
	if ((vm.count("stepwidth-algorithm"))
			&& (algorithm_name != MinimizerFactory::getNameForType(MinimizerFactory::landweber))) {
		LOG(error, "stepwidth algorithm set but not Landweber chosen.");
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_stopping_criteria() const
{
	if ((vm.count("stopping-criteria")) && (stopping_criteria.empty())) {
		LOG(error, "Empty stopping criteria set.");
		return false;
	}
	return true;
}

bool CommandLineOptions::checkSensibility_updatealgorithm() const
{
	if (vm.count("update-algorithm"))
		if (!((updatetype >= LastNSearchDirections::RoundRobin)
				&& (updatetype < LastNSearchDirections::MAX_UpdateAlgorithmType))) {
			LOG(error, "Illegal update type set.");
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
			LOG(error, "Illegal Wolfe constants given, must be two and positive.");
			return false;
		}
	return true;
}

bool CommandLineOptions::checkSensibility() const
{
	Options::checkSensibility();

	bool status = true;
	status &= checkSensibility_config();
	status &= checkSensibility_delta();
	status &= checkSensibility_everynthtuple();
	status &= checkSensibility_OrthogonalDirections();
	status &= checkSensibility_regularizationparameter();
	status &= checkSensibility_tau();
	status &= checkSensibility_tuple_parameters();
	status &= checkSensibility_algorithm();
	status &= checkSensibility_minlib();
	status &= checkSensibility_norms();
	status &= checkSensibility_pvalues();
	status &= checkSensibility_searchspace();
	status &= checkSensibility_stepwidth_algorithm();
	status &= checkSensibility_stopping_criteria();
	status &= checkSensibility_updatealgorithm();
	status &= checkSensibility_wolfeconstants();
	status &= internal_checkSensibility();

	if (!status)
		showHelpinErrorCase();

	return status;
}

void CommandLineOptions::setSecondaryValues()
{
	Options::setSecondaryValues();

	// set search space related stuff
	SearchspaceFactory::setCurrentType(
			SearchspaceFactory::getType(searchspace_type)
	);
	FunctionalMinimizerFactory::setMinLib(minlib);

	type = MinimizerFactory::getTypeForName(algorithm_name);

	// set good maximum of inner iterations
	if (maxinneriter == 0)
		maxinneriter = N*100;

	if (stopping_criteria.empty()) {
		int count = 0;
		if (vm.count("delta")) {
			if (count != 0)
				stopping_criteria += " || ";
			if (vm.count("algorithm") && (algorithm_name == "RESESOP"))
				stopping_criteria += "Residuum";
			else
				stopping_criteria += "RelativeResiduum";
			++count;
		}
		if (vm.count("maxiter")) {
			if (count != 0)
				stopping_criteria += " || ";
			stopping_criteria += "MaxIterationCount";
			++count;
		}
		if (vm.count("max-walltime")) {
			if (count != 0)
				stopping_criteria += " || ";
			stopping_criteria += "MaxWalltime";
			++count;
		}
	}

	// set lp value in case of regularized_l1
	if ((type_spacex == "regularized_l1") && (vm.count("px") == 0))
		px = 1.;
	if ((type_spacey == "regularized_l1") && (vm.count("py") == 0))
		py = 1.;

	internal_setSecondaryValues();
}

void CommandLineOptions::store(std::ostream &_output) const
{
	Options::store(_output);

	_output << "# [General]" << std::endl;
	writeValue<bool>(_output, vm,  "database-replace");
	writeValue<unsigned int>(_output, vm,  "output-steps");
	if (tuple_parameters.size() != 0) {
		for (std::vector<std::string>::const_iterator iter = tuple_parameters.begin();
				iter != tuple_parameters.end();)
			_output << "\ttuple-parameters = " << *(iter++) << std::endl;
	}

	_output << "# [Banach Space]" << std::endl;
	writeValue<std::string>(_output, vm,  "type-space-x");
	writeValue<std::string>(_output, vm,  "type-space-y");
	writeValue<double>(_output, vm,  "px");
	writeValue<double>(_output, vm,  "py");
	writeValue<double>(_output, vm,  "powerx");
	writeValue<double>(_output, vm,  "powery");
	writeValue<double>(_output, vm,  "regularization-parameter");

	_output << "# [Algorithm]" << std::endl;
	writeValue<std::string>(_output, vm,  "algorithm");
	writeValue<double>(_output, vm,  "delta");
	writeValue<unsigned int>(_output, vm,  "every-nth-tuple");
	writeValue<boost::filesystem::path>(_output, vm, "iteration-file");
	writeValue<std::string>(_output, vm,  "minimization-library");
	writeValue<std::string>(_output, vm,  "stopping-criteria");
	writeValue<unsigned int>(_output, vm,  "max-inner-iterations");
	writeValue<unsigned int>(_output, vm,  "maxiter");
	writeValue<unsigned int>(_output, vm,  "max-walltime");
	writeValue<double>(_output, vm,  "tau");

	_output << "# [Landweber]" << std::endl;
	writeValue<double>(_output, vm,  "C");
	writeValue<unsigned int>(_output, vm,  "stepwidth-algorithm");

	_output << "# [SESOP]" << std::endl;
	writeValue<unsigned int>(_output, vm,  "number-directions");
	writeValue<bool>(_output, vm,  "calculateAngles");
	writeValue<bool>(_output, vm,  "enforceRandomMapping");
	writeValue<bool>(_output, vm,  "inexact-linesearch");
	writeValue<unsigned int>(_output, vm,  "orthogonal-directions");
	writeValue<std::string>(_output, vm,  "searchspace");
	writeValue<unsigned int>(_output, vm,  "update-algorithm");
	if (vm.count("wolfe-constants") != 0) {
		const std::vector<double> values =
				vm["wolfe-constants"].as< std::vector<double> >();
		for (std::vector<double>::const_iterator iter = values.begin();
				iter != values.end();)
			_output << "\twolfe-constants = " << *(iter++) << std::endl;
	}

	_output << "# [Split Feasibility Problem]" << std::endl;
	writeValue<std::string>(_output, vm,  "auxiliary-constraints");
	writeValue<unsigned int>(_output, vm,  "max-sfp-loops");

	internal_store(_output);
}
