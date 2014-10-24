/*
 * SmoothnessModulusUnitTest.cpp
 *
 *  Created on: Apr 03, 2014
 *      Author: heber
 */

#include "SmoothnessModulusUnitTest.hpp"

#include <cassert>
#include <Eigen/Dense>

#include "Minimizations/MinimizationExceptions.hpp"
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
