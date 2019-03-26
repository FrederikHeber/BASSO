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

    class_<SpaceElement, SpaceElement_ptr_t >(
    		"SpaceElement",
			"Element of a NormedSpace",
			no_init)
        .def("isInSpace", &SpaceElement::isInSpace,
        		"Predicate for whether the element is in the given otherspace or not.",
				args("self", "otherspace"))
        .def("isZero", SpaceElement_isZero,
        		"Predicate whether this element has only zero components",
				args("self"))
        .def("isZero", SpaceElement_isZero_tolerance,
        		"Predicate whether this element has only approximately zero components up to tolerance in magnitude",
        		args("self", "tolerance"))
        .def("isApproxToConstant", &SpaceElement::isApproxToConstant,
        		"Predicate whether this element has only approximately constant components up to tolerance in magnitude",
        		args("self", "constant", "tolerance"))
        .def("isApprox", &SpaceElement_isApprox,
        		"Predicate whether this element is approximately equivalent to otherelement",
				args("self", "otherelement", "tolerance"))
        .def("isNonnegative", &SpaceElement::isNonnegative,
        		"Predicate whether this element has only non negative components",
				args("self"))
        .def("getSignVector", &SpaceElement::getSignVector,
        		"Returns a new vector with each component as the sign of this element's component, i.e. -1,0,1",
				args("self"))
        .def("getAbsVector", &SpaceElement::getAbsVector,
        		"Returns a new vector with each component as the absolute value",
				args("self"))
        .def("getCircShiftedVector", &SpaceElement::getCircShiftedVector,
        		"Returns a new vector where components have been shifted circularly by shift entries",
				args("self", "shift"))
        .def("getMaxCoefficientAndIndex", &SpaceElement::getMaxCoefficientAndIndex,
        		"Getter for the pair of the maximal component and its index",
				args("self"))
        .def_readonly("norm", &SpaceElement::Norm,
        		"Returns the norm of this element")
        .def("pow", &SpaceElement::pow,
        		"Calculates each element to the power of the given exponent",
				args("self", "exponent"))
        .def_readonly("space", &SpaceElement::getSpace,
        		"Getter for the NormedSpace instance this element is associated to")
        .def("setZero", &SpaceElement::setZero,
        		"Sets all components to zero",
				args("self"))
        .def("__getitem__", &SpaceElement_getitem,
        		"Getter for a single component with index",
				args("self", "index"))
        .def("__setitem__", &SpaceElement_setitem,
        		"Setter for a single component with index to value",
				args("self", "index", "value"))
        .def("get", &pyBasso_SpaceElement_access::get,
        		"Returns the SpaceElement as a minieigen.VectorX",
				args("self"))
        .def("set", &pyBasso_SpaceElement_access::set,
        		"Setter of this SpaceElement from the components in vector of type minieigen.VectorX",
				args("self", "vector"))
        .def(self_ns::str(self_ns::self))
        .def(self_ns::self * float())
        .def(self_ns::self + SpaceElement_ptr_t())
        .def(self_ns::self += SpaceElement_ptr_t())
        .def(self_ns::self - SpaceElement_ptr_t())
        .def(self_ns::self -= SpaceElement_ptr_t())
        .def(self_ns::self == self_ns::self)
    ;
}
