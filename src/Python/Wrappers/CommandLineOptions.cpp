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

void export_commandlineoptions()
{
    class_<CommandLineOptions>("CommandLineOptions", init<>())
        .def_readwrite("algorithm_name", &CommandLineOptions::algorithm_name)
        .def_readwrite("auxiliary_constraints", &CommandLineOptions::auxiliary_constraints)
        .def_readwrite("C", &CommandLineOptions::C)
        .def_readwrite("calculateAngles", &CommandLineOptions::calculateAngles)
        .def_readwrite("delta", &CommandLineOptions::delta)
        .def_readwrite("enforceRandomMapping", &CommandLineOptions::enforceRandomMapping)
        .def_readwrite("everynthtuple", &CommandLineOptions::everynthtuple)
        .def_readwrite("inexactLinesearch", &CommandLineOptions::inexactLinesearch)
        .def_readwrite("iteration_file", &CommandLineOptions::iteration_file)
        .def_readwrite("maxinneriter", &CommandLineOptions::maxinneriter)
        .def_readwrite("maxiter", &CommandLineOptions::maxiter)
        .def_readwrite("max_sfp_loops", &CommandLineOptions::max_sfp_loops)
        .def_readwrite("maxwalltime", &CommandLineOptions::maxwalltime)
        .def_readwrite("minlib", &CommandLineOptions::minlib)
        .def_readwrite("N", &CommandLineOptions::N)
        .add_property("orthogonalization_type",
                      &CommandLineOptions_orthogonalization_type_get, // getter
                      &CommandLineOptions_orthogonalization_type_set) // setter
        .def_readwrite("outputsteps", &CommandLineOptions::outputsteps)
        .def_readwrite("regularization_parameter", &CommandLineOptions::regularization_parameter)
        .def_readwrite("searchspace_type", &CommandLineOptions::searchspace_type)
        .def_readwrite("stepwidth_type", &CommandLineOptions::stepwidth_type)
        .def_readwrite("stopping_criteria", &CommandLineOptions::stopping_criteria)
        .def_readwrite("tau", &CommandLineOptions::tau)
        .def_readwrite("tolerance_linesearch", &CommandLineOptions::tolerance_linesearch)
        .def_readwrite("tolerance_spacex", &CommandLineOptions::tolerance_spacex)
        .def_readwrite("tuple_parameters", &CommandLineOptions::tuple_parameters)
        .def_readwrite("updatetype", &CommandLineOptions::updatetype)
        .def_readwrite("wolfe_constants", &CommandLineOptions::wolfe_constants)
        .def("__repr__", &CommandLineOptions_toString)
        .def("setValues", &CommandLineOptions::setSecondaryValues)
        .def("checkSensibility", &CommandLineOptions::checkSensibility)
        .def("setVerbosity", &setVerbosity)
        ;
}

