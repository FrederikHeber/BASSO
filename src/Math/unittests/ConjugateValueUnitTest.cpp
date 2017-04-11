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
 * ConjugateValueUnitTest.cpp
 *
 *  Created on: Oct 28, 2014
 *      Author: heber
 */

#include "ConjugateValueUnitTest.hpp"

#include <cmath>
#include <limits>

#include "Math/Helpers.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( ConjugateValueUnitTest );


void ConjugateValueUnitTest::setUp()
{
}


void ConjugateValueUnitTest::tearDown()
{
}

void ConjugateValueUnitTest::Test()
{
	CPPUNIT_ASSERT( fabs(2. - Helpers::ConjugateValue(2.)) < std::numeric_limits<double>::epsilon()*1e2);
	CPPUNIT_ASSERT( fabs(1.5 - Helpers::ConjugateValue(3.)) < std::numeric_limits<double>::epsilon()*1e2 );
	CPPUNIT_ASSERT( fabs(3. - Helpers::ConjugateValue(1.5)) < std::numeric_limits<double>::epsilon()*1e2 );
	CPPUNIT_ASSERT( fabs(4. - Helpers::ConjugateValue(4./3.)) < std::numeric_limits<double>::epsilon()*1e2 );
	CPPUNIT_ASSERT( fabs((10./9.) - Helpers::ConjugateValue(10.)) < std::numeric_limits<double>::epsilon()*1e2 );
}
