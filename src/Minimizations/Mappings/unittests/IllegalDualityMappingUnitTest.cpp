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
 * IllegalDualityMappingUnitTest.cpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "IllegalDualityMappingUnitTest.hpp"

#include <boost/assign.hpp>

#include <Eigen/Dense>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/IllegalDualityMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( IllegalDualityMappingUnitTest );

using namespace boost::assign;

void IllegalDualityMappingUnitTest::setUp()
{
}


void IllegalDualityMappingUnitTest::tearDown()
{
}

void IllegalDualityMappingUnitTest::IllegalCall()
{
	IllegalDualityMapping J_illegal;
	Eigen::VectorXd xtemp(10);
	xtemp << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929;
	NormedSpaceFactory::args_t args;
	args += boost::any(2.), boost::any(2.);
	const NormedSpace_ptr_t SpaceX =
			NormedSpaceFactory::create(2, "lp", args);
	SpaceElement_ptr_t x = SpaceX->createElement();
	CPPUNIT_ASSERT_THROW( J_illegal(x), MinimizationIllegalValue_exception);
}
