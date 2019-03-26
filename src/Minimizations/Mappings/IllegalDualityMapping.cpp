/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
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
 * IllegalDualityMapping.cpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#include "IllegalDualityMapping.hpp"

#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

IllegalDualityMapping::IllegalDualityMapping() :
	LpDualityMapping(NormedSpaceFactory::getDummySpace(), 2.)
{}

const NormedSpace_ptr_t IllegalDualityMapping::getSourceSpace() const
{
	// we just throw as this function must now be called
	throw MinimizationIllegalValue_exception()
		<< MinimizationIllegalValue_name("illegally called");
}

const NormedSpace_ptr_t IllegalDualityMapping::getTargetSpace() const
{
	// we just throw as this function must now be called
	throw MinimizationIllegalValue_exception()
		<< MinimizationIllegalValue_name("illegally called");
}

void IllegalDualityMapping::operator()(
		const SpaceElement_ptr_t &_sourceelement,
		SpaceElement_ptr_t &_destelement
		) const
{
	// we just throw as this function must now be called
	throw MinimizationIllegalValue_exception()
		<< MinimizationIllegalValue_name("illegally called");
}

const Mapping_ptr_t IllegalDualityMapping::getAdjointMapping() const
{
	// we just throw as this function must now be called
	throw MinimizationIllegalValue_exception()
		<< MinimizationIllegalValue_name("illegally called");
}


