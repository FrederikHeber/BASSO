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
 * SpaceElement.cpp
 *
 *  Created on: Aug 9, 2018
 *      Author: heber
 */

#include <boost/python.hpp>
#include <boost/python/operators.hpp>

#include <string>

#include "Minimizations/Elements/SpaceElement.hpp"

#include "Python/utility.hpp"

using namespace boost::python;

void export_spaceelement()
{
    const bool (SpaceElement::*SpaceElement_isZero_tolerance)(const double) const =
            &SpaceElement::isZero;
    SpaceElement_ptr_t (*double_times_SpaceElement)(
            const double _alpha,
            const SpaceElement_ptr_t &_element) = &operator*;
    class_<SpaceElement, SpaceElement_ptr_t >("SpaceElement", no_init)
        .def("isInSpace", &SpaceElement::isInSpace)
        .def("isZero", SpaceElement_isZero)
        .def("isZero", SpaceElement_isZero_tolerance)
        .def("isApproxToConstant", &SpaceElement::isApproxToConstant)
        .def("isApprox", &SpaceElement_isApprox)
        .def("isNonnegative", &SpaceElement::isNonnegative)
        .def("getSignVector", &SpaceElement::getSignVector)
        .def("getAbsVector", &SpaceElement::getAbsVector)
        .def("getCircShiftedVector", &SpaceElement::getCircShiftedVector)
        .def("getMaxCoefficientAndIndex", &SpaceElement::getMaxCoefficientAndIndex)
        .def_readonly("norm", &SpaceElement::Norm)
        .def("pow", &SpaceElement::pow)
        .def_readonly("space", &SpaceElement::getSpace)
        .def("setZero", &SpaceElement::setZero)
        .def("__getitem__", &SpaceElement_getitem)
        .def("__setitem__", &SpaceElement_setitem)
        .def("get", &pyBasso_SpaceElement_access::get)
        .def("set", &pyBasso_SpaceElement_access::set)
        .def(self_ns::str(self_ns::self))
        .def(self_ns::self * float())
        .def(self_ns::self + SpaceElement_ptr_t())
        .def(self_ns::self += SpaceElement_ptr_t())
        .def(self_ns::self - SpaceElement_ptr_t())
        .def(self_ns::self -= SpaceElement_ptr_t())
        .def(self_ns::self == self_ns::self)
    ;
}
