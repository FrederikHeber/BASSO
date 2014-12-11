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
