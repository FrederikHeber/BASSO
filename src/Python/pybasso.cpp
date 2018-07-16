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

/* pybasso.cpp
 *
 *  Created on: Jul 16, 2018
 *      Author: heber
 */

#include <boost/assign.hpp>
#include <boost/python.hpp>

#include <string>

#include "Minimizations/InverseProblems/InverseProblemFactory.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

using namespace boost::assign;
using namespace boost::python;

NormedSpace_ptr_t create_LpSpace(
		const int _dimension,
		const double _p,
		const double _power
		)
{
	InverseProblemFactory::args_t _args_SpaceX;
	_args_SpaceX += boost::any(_p), boost::any(_power);
	return NormedSpaceFactory::create(
			_dimension, std::string("lp"), _args_SpaceX);
}

BOOST_PYTHON_MODULE(pyBasso)
{
	// classes
	class_<LpNorm>("LpNorm", no_init)
		.def("getPvalue", &LpNorm::getPvalue)
	;

    class_<NormedSpace>("NormedSpace", no_init)
        .def("getDimension", &NormedSpace::getDimension)
    ;

    // factory methods
    def("create_NormedSpace", &NormedSpaceFactory::create);
}
