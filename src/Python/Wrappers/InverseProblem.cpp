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
 * InverseProblem.cpp
 *
 *  Created on: Aug 9, 2018
 *      Author: heber
 */

#include <boost/python.hpp>

#include "Minimizations/InverseProblems/InverseProblem.hpp"

#include "Python/utility.hpp"

using namespace boost::python;

void export_inverseproblem()
{
    class_<InverseProblem, InverseProblem_ptr_t>(
    		"InverseProblem",
			"Inverse problem instance that contains all objects for solving problems of the form F(x)=y",
			init<
                const Mapping_ptr_t, const NormedSpace_ptr_t,
                const NormedSpace_ptr_t, const SpaceElement_ptr_t>())
        .def_readonly("sourcespace", &InverseProblem::SourceSpace,
        		"Source NormedSpace of the mapping A")
        .def_readonly("targetspace", &InverseProblem::TargetSpace,
        		"Target NormedSpace of the mapping A")
        .def_readonly("dualsourcespace", &InverseProblem::DualSourceSpace,
        		"Dual space of the source space")
        .def_readonly("dualtargetspace", &InverseProblem::DualTargetSpace,
        		"Dual space of the target space")
        .def_readonly("A", &InverseProblem::A,
        		"Linear or non-linear Mapping to invert")
        .def_readonly("At", &InverseProblem::A_t,
        		"Adjoint of the mapping A")
        .def_readwrite("x", &InverseProblem::x,
        		"Solution of the inverse problem")
        .def_readonly("y", &InverseProblem::y,
        		"Right-hand side y of the inverse problem")
        .def("solve", &InverseProblem_solve,
        		"Function to solve the inverse problem")
        .def("project", &RangeProjection_solve,
        		"Function to project the right-hand side onto the range of  the mapping")
    ;
}


