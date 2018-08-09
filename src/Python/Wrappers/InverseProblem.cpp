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
    class_<InverseProblem, InverseProblem_ptr_t>("InverseProblem", init<
                const Mapping_ptr_t, const NormedSpace_ptr_t,
                const NormedSpace_ptr_t, const SpaceElement_ptr_t>())
        .def_readonly("sourcespace", &InverseProblem::SourceSpace)
        .def_readonly("targetspace", &InverseProblem::TargetSpace)
        .def_readonly("dualsourcespace", &InverseProblem::DualSourceSpace)
        .def_readonly("dualtargetspace", &InverseProblem::DualTargetSpace)
        .def_readonly("A", &InverseProblem::A)
        .def_readonly("At", &InverseProblem::A_t)
        .def_readwrite("x", &InverseProblem::x)
        .def_readonly("y", &InverseProblem::y)
        .def("solve", &InverseProblem_solve)
        .def("project", &RangeProjection_solve)
    ;
}


