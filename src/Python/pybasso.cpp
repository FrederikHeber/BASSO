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

#include <boost/python.hpp>

void export_commandlineoptions();
void export_inverseproblem();
void export_mapping();
void export_norm();
void export_normedspace();
void export_spaceelement();

BOOST_PYTHON_MODULE(pyBasso)
{
    /*** classes with virtual functions ***/
	export_norm();
	export_mapping();

    /*** classes ***/
	export_normedspace();
	export_spaceelement();
	export_inverseproblem();
	export_commandlineoptions();
}
