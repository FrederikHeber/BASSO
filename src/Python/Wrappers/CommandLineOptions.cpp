/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2018 Frederik Heber. All rights reserved.
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
 * CommandLineOptions.cpp
 *
 *  Created on: Aug 9, 2018
 *      Author: heber
 */

#include <boost/python.hpp>

#include "Options/CommandLineOptions.hpp"

#include "Python/utility.hpp"

using namespace boost::python;

void setVerbosity(CommandLineOptions &opts, unsigned int verbose)
{
	opts.verbose = verbose;
	opts.setVerbosity();
}

void CommandLineOptions_tuple_parameters_set(CommandLineOptions &opts, list &_list)
{
	for (int i = 0; i < len(_list); ++i)
		opts.tuple_parameters.push_back(extract<std::string>(_list[i]));
}

boost::python::list CommandLineOptions_tuple_parameters_get(CommandLineOptions &opts)
{
	boost::python::list result;
	// call original function
	std::vector<std::string> v = opts.tuple_parameters;
	// put all the strings inside the python list
	std::vector<std::string>::iterator it;
	for (it = v.begin(); it != v.end(); ++it){
	   result.append(*it);
	}
	return result;
}


void export_commandlineoptions()
{
    class_<CommandLineOptions>(
    		"Options",
			"Options structure that control the optimization",
			init<>(args("self")))
        .def_readwrite("algorithm_name", &CommandLineOptions::algorithm_name,
        		"Name of optimization algorithm to use: Landweber, SESOP, RESESOP")
        .def_readwrite("auxiliary_constraints", &CommandLineOptions::auxiliary_constraints,
        		"Auxiliary constraints (Nonnegative, Nonpositive, Unity) with boolean AND (&&)")
        .def_readwrite("C", &CommandLineOptions::C,
        		"C value")
        .def_readwrite("calculateAngles", &CommandLineOptions::calculateAngles,
        		"Whether the calculate angles between search directions")
        .def_readwrite("delta", &CommandLineOptions::delta,
        		"delta threshold for stopping criterion")
        .def_readwrite("enforceRandomMapping", &CommandLineOptions::enforceRandomMapping,
        		"whether to enforce random mapping or not")
        .def_readwrite("everynthtuple", &CommandLineOptions::everynthtuple,
        		"")
        .def_readwrite("inexactLinesearch", &CommandLineOptions::inexactLinesearch,
        		"Whether to use Wolfe criterion for performing an inexact line search")
        .def_readwrite("iteration_file", &CommandLineOptions::iteration_file,
        		"SQLite file containing information on the iteration after optimization")
        .def_readwrite("maxinneriter", &CommandLineOptions::maxinneriter,
        		"Maximum number of inner iterations in split feasibility problems")
        .def_readwrite("maxiter", &CommandLineOptions::maxiter,
        		"Maximum number of iterations")
        .def_readwrite("max_sfp_loops", &CommandLineOptions::max_sfp_loops,
        		"Maximum number of Split Feasibility Loops iteration")
        .def_readwrite("maxwalltime", &CommandLineOptions::maxwalltime,
        		"Maximum walltime to use for optimization as stopping criterion")
        .def_readwrite("minlib", &CommandLineOptions::minlib,
        		"Name of the minimization library: gsl, nlopt")
        .def_readwrite("N", &CommandLineOptions::N,
        		"Number of search directions")
        .add_property("orthogonalization_type",
                      &CommandLineOptions_orthogonalization_type_get, // getter
                      &CommandLineOptions_orthogonalization_type_set) // setter
        .def_readwrite("outputsteps", &CommandLineOptions::outputsteps,
        		"")
        .def_readwrite("regularization_parameter", &CommandLineOptions::regularization_parameter,
        		"Parameter for regularization")
        .def_readwrite("searchspace_type", &CommandLineOptions::searchspace_type,
        		"Search space type to use: LastNDirections, Nemirovsky")
        .def_readwrite("stepwidth_type", &CommandLineOptions::stepwidth_type,
        		"Type of step width choosing procedure: LandweberFixes(0), MinimizingResidual(1) using Brent's method, ConstantRegularizedL1Norm (2), DynamicRegularizedL1Norm (3)")
        .def_readwrite("stopping_criteria", &CommandLineOptions::stopping_criteria,
        		"Stopping criterion (DivergentResiduum, MaxIterationCount, RelativeChangeResiduum, RelativeResiduum, Residuum, MaxWalltime) with boolean combination,\ne.g. RelativeResiduum || MaxWalltime")
        .def_readwrite("tau", &CommandLineOptions::tau,
        		"Tolerance parameter in case of noisy data")
        .def_readwrite("tolerance_linesearch", &CommandLineOptions::tolerance_linesearch,
        		"Tolerance threshold to use for stopping linesearch")
        .def_readwrite("tolerance_spacex", &CommandLineOptions::tolerance_spacex,
        		"Tolerance threshold to use for duality mapping between spaces")
        .add_property("tuple_parameters",
        		&CommandLineOptions::tuple_parameters,
				&CommandLineOptions_tuple_parameters_set,
        		"Tuple of parameters to add to the table written to the iteration_file")
        .def_readwrite("updatetype", &CommandLineOptions::updatetype,
        		"How to iterate through the set of search directions: RoundRobin(0), MostParallel(1), MostOrthogonal(2)")
        .def_readwrite("wolfe_constants", &CommandLineOptions::wolfe_constants,
        		"")
        .def("__repr__", &CommandLineOptions_toString,
        		"String representation of the Options, NOT WORKING at the moment")
        .def("setValues", &CommandLineOptions::setSecondaryValues,
        		"Function to set all secondary values, needs to be called after changes",
				args("self"))
        .def("checkSensibility", &CommandLineOptions::checkSensibility,
        		"Function to check sensibility of optimization values",
				args("self"))
        .def("setVerbosity", &setVerbosity,
        		"Setting the verbosity of the optimization",
				args("level"))
        ;
}

