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
    class_<NormedSpace, NormedSpace_ptr_t >(
    		"NormedSpace",
			"Normed space instance that is inherently connected with a Norm instance",
			no_init)
    /* does not work: .def_readonly("norm", &NormedSpace_getNorm) */
        .def("getNorm", &NormedSpace_getNorm, return_internal_reference<1>(),
        		"getter for the internal norm object")
        .def_readonly("space", &NormedSpace::getSpace,
        		"getter for the space instance")
        .def_readonly("dualspace", &NormedSpace::getDualSpace,
        		"getter for the dual space")
        .def("getDualityMapping", &NormedSpace::getDualityMapping,
        		"getter for the internal duality mapping between this space and its dual")
        .def_readonly("dim", &NormedSpace::getDimension,
        		"Getter for the finite dimension of this vector space")
        .def("createElement", &NormedSpace::createElement,
        		"Creates a SpaceElement instance associated with this NormedSpace")
    ;

    def("create_LpSpace", &create_LpSpace,
    		"Creator for an Lp space instance with an Lp Norm");
}

