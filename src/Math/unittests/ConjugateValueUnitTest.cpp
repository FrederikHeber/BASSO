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
