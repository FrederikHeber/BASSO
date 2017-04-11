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
 * NormFactoryUnitTest.cpp
 *
 *  Created on: Dec 11, 2014
 *      Author: heber
 */

#include "NormFactoryUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/Norms/NormExceptions.hpp"
#include "Minimizations/Norms/NormFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( NormFactoryUnitTest );


void NormFactoryUnitTest::setUp()
{
}


void NormFactoryUnitTest::tearDown()
{
}

void NormFactoryUnitTest::throwTest()
{
	// we check that assertion is thrown for invalid p value
//	std::cout << "The following assertion is intended and does not indicate a failure of the test." << std::endl;
	CPPUNIT_ASSERT_THROW(
			NormFactory::createLpInstance(-0.5),
			NormIllegalValue_exception );
}
