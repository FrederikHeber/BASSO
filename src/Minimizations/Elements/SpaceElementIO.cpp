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
 * SpaceElementIO.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SpaceElementIO.hpp"

#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

void SpaceElementWriter::output(std::ostream &_ost, const SpaceElement_ptr_t &_element)
{
	using namespace MatrixIO;
	_ost << RepresentationAdvocate::get(_element);
}

void SpaceElementReader::input(std::istream &_ist, SpaceElement_ptr_t &_element)
{
	using namespace MatrixIO;
	MatrixIO::dims_t dims = getDimensions(_ist);
	// assert its a row vector
	assert( _element->getSpace()->getDimension() == dims.first );
	assert( 1 == dims.second );
	Eigen::VectorXd vector;
	_ist >> vector;
	RepresentationAdvocate::set(_element, vector);
}
