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
 * Norm.cpp
 *
 *  Created on: Aug 9, 2018
 *      Author: heber
 */

#include <boost/python.hpp>

#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"

using namespace boost::python;

struct NormWrap : Norm, wrapper<Norm> {
    virtual bool isSmooth() const
    { return this->get_override("isSmooth")(); }
    virtual const double internal_operator(const SpaceElement_ptr_t &_element) const
    { return this->get_override("internal_operator")(); }
};

void export_norm()
{
    class_<NormWrap, Norm_ptr_t, boost::noncopyable>(
    		"Norm",
    		"Norm instance to measure distances, inherently bound to a NormedSpace.\nUse create_LpSpace() to create Norm and Space together.",
			no_init)
        .def_readonly("p", &Norm::getPvalue,
        		"p value of the norm")
        .def_readonly("space", &Norm::getSpace,
        		"reference to the associated NormedSpace of this norm")
        .def("__call__", &Norm::operator(),
        		"Calculating the norm of the SpaceElement x",
				args("self", "x"))
        ;
    class_<LpNorm, bases<Norm> >(
    		"LpNorm",
			"Lp norm for spaces of sequences",
			no_init);
}
