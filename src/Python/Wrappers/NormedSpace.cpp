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
 * NormedSpace.cpp
 *
 *  Created on: Aug 9, 2018
 *      Author: heber
 */

#include <boost/python.hpp>

#include "Minimizations/Spaces/NormedSpace.hpp"

#include "Python/utility.hpp"

using namespace boost::python;

void export_normedspace()
{
    class_<NormedSpace, NormedSpace_ptr_t >("NormedSpace", no_init)
    /* does not work: .def_readonly("norm", &NormedSpace_getNorm) */
        .def("getNorm", &NormedSpace_getNorm, return_internal_reference<1>())
        .def_readonly("space", &NormedSpace::getSpace)
        .def_readonly("dualspace", &NormedSpace::getDualSpace)
        .def("getDualityMapping", &NormedSpace::getDualityMapping)
        .def_readonly("dim", &NormedSpace::getDimension)
        .def("createElement", &NormedSpace::createElement)
    ;

    def("create_LpSpace", &create_LpSpace);
}

