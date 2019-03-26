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
 * SmoothnessModulusUnitTest.cpp
 *
 *  Created on: Apr 03, 2014
 *      Author: heber
 */

#include "SmoothnessModulusUnitTest.hpp"

#include <cassert>
#include <Eigen/Dense>

#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Functions/SmoothnessModulus.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SmoothnessModulusUnitTest );


void SmoothnessModulusUnitTest::setUp()
{
}


void SmoothnessModulusUnitTest::tearDown()
{
}

void SmoothnessModulusUnitTest::throwTest()
{
	// we check that assertion is thrown for invalid p value
//	std::cout << "The following assertion is intended and does not indicate a failure of the test." << std::endl;
	CPPUNIT_ASSERT_THROW(
			SmoothnessModulus modulus(0.5),
			MinimizationIllegalValue_exception );
}

void SmoothnessModulusUnitTest::oneandhalftest()
{
	SmoothnessModulus modulus(1.5);

//	CPPUNIT_ASSERT_EQUAL( 0.235702260395516, modulus( 0.5));
	CPPUNIT_ASSERT( fabs(0.235702260395516 - modulus( 0.5)) < 1e-10 );
}

void SmoothnessModulusUnitTest::twotest()
{
	SmoothnessModulus modulus(2.);

//	CPPUNIT_ASSERT_EQUAL( 0.125, modulus( 0.5));
	CPPUNIT_ASSERT( fabs(0.125 - modulus( 0.5)) < 1e-10 );
}

void SmoothnessModulusUnitTest::threetest()
{
	SmoothnessModulus modulus(3.);

//	CPPUNIT_ASSERT_EQUAL( 0.25, modulus( 0.5));
	CPPUNIT_ASSERT( fabs(0.25 - modulus( 0.5)) < 1e-10 );
}
