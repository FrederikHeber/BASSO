/*
 * BregmanDistanceUnitTest.cpp
 *
 *  Created on: Apr 25, 2014
 *      Author: heber
 */

#include "BregmanDistanceUnitTest.hpp"

#include <Eigen/Dense>

#include "Log/Logging.hpp"
#include "Minimizations/BregmanDistance.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( BregmanDistanceUnitTest );


void BregmanDistanceUnitTest::setUp()
{
	// BregmanDistance uses logging
	startLogging();
}


void BregmanDistanceUnitTest::tearDown()
{
}

void BregmanDistanceUnitTest::oneNorm()
{
	CPPUNIT_ASSERT_THROW( BregmanDistance d_1(1),
			MinimizationIllegalValue_exception);
	Eigen::VectorXd t(2);
	t << 3.5,5;
	Eigen::VectorXd x(2);
	x << 4,3;
//	CPPUNIT_ASSERT_EQUAL( -52.5, d_1(t,x) );
//	std::cout << "BregmanDistance d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}

void BregmanDistanceUnitTest::twoNorm()
{
	BregmanDistance d_2(2);
	Eigen::VectorXd t(2);
	t << 3.5,5;
	Eigen::VectorXd x(2);
	x << 4,3;
//	CPPUNIT_ASSERT_EQUAL( 2.125, d_2(t,x,q) );
	CPPUNIT_ASSERT( fabs( 2.125 - d_2(t,x)) < BASSOTOLERANCE);
//	std::cout << "BregmanDistance d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}

void BregmanDistanceUnitTest::inftyNorm()
{
	CPPUNIT_ASSERT_THROW( BregmanDistance d_1(1),
			MinimizationIllegalValue_exception);
	Eigen::VectorXd t(2);
	t << 3.5,5;
	Eigen::VectorXd x(2);
	x << 4,3;
//	CPPUNIT_ASSERT_EQUAL( -6, d_infty(t,x) );
//	std::cout << "BregmanDistance d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}
