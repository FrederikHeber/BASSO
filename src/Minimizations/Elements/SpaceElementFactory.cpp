/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
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
 * SpaceElementFactory.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

#include "MatrixIO/MatrixIO.hpp"
#include "MatrixIO/MatrixIOExceptions.hpp"
#include "Minimizations/Elements/SpaceElementFactory.hpp"
#include "Minimizations/Elements/SpaceElementIO.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

SpaceElement_ptr_t SpaceElementFactory::create(
		const NormedSpace_weakptr_t _SpaceRef,
		const std::string &_name
		)
{
	SpaceElement_ptr_t vector =
			NormedSpace_ptr_t(_SpaceRef)->createElement();

	using namespace MatrixIO;
	if (MatrixIO::isPresentFile(_name)) {
		std::ifstream ist(_name.c_str());
		if (ist.good()) {
			try {
				SpaceElementReader::input(ist, vector);
			} catch (MatrixIOStreamEnded_exception &e) {
				std::cerr << "Failed to fully parse vector from " << _name << std::endl;
			}
		} else {
			std::cerr << "Failed to open " << _name << std::endl;
		}
	} else {
		std::cerr << _name << " is not a valid name for SpaceElementFactory." << std::endl;
	}

	return vector;
}
