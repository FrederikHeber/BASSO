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
 * NormedSpace.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormedSpace.hpp"

#include <cassert>
#include <boost/bind.hpp>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/VectorSpaceOperationCounts.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

NormedSpace::NormedSpace(
		const unsigned int _dimension) :
	dimension(_dimension)
{}

NormedSpace::NormedSpace(
		const unsigned int _dimension,
		const Norm_ptr_t &_norm
		) :
	norm(_norm),
	dimension(_dimension)
{}

SpaceElement_ptr_t NormedSpace::createElement() const
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::ElementCreation>(
					opcounts.instance));
	SpaceElement_ptr_t newelement(new SpaceElement(Space));
	newelement->setSelfRef(newelement);
	return newelement;
}

const Mapping_ptr_t& NormedSpace::getDualityMapping() const
{
	assert (dualitymapping != NULL);
	return dualitymapping;
}
