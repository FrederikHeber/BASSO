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
	LpNorm lpnorm(1.);
	DualityMapping J_1(1.);
	CPPUNIT_ASSERT_THROW( BregmanDistance d_1(lpnorm, J_1, 1),
			MinimizationIllegalValue_exception);
	Eigen::VectorXd t(2);
	t << 3.5,5;
	Eigen::VectorXd x(2);
	x << 4,3;
//	CPPUNIT_ASSERT_EQUAL( -52.5, d_1(t,x, 1) );
//	std::cout << "BregmanDistance d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}

void BregmanDistanceUnitTest::twoNorm()
{
	LpNorm lpnorm(2.);
	DualityMapping J_2(2.);
	BregmanDistance d_2(lpnorm, J_2, 2);
	Eigen::VectorXd t(2);
	t << 3.5,5;
	Eigen::VectorXd x(2);
	x << 4,3;
//	CPPUNIT_ASSERT_EQUAL( 2.125, d_2(t,x,q) );
	CPPUNIT_ASSERT( fabs( 2.125 - d_2(t,x,2.)) < BASSOTOLERANCE);
//	std::cout << "BregmanDistance d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}

void BregmanDistanceUnitTest::inftyNorm()
{
	LpNorm lpnorm(LpNorm::Infinity);
	DualityMapping J_infty(LpNorm::Infinity);
	CPPUNIT_ASSERT_THROW( BregmanDistance d_infty(lpnorm, J_infty, 0),
			MinimizationIllegalValue_exception);
	Eigen::VectorXd t(2);
	t << 3.5,5;
	Eigen::VectorXd x(2);
	x << 4,3;
//	CPPUNIT_ASSERT_EQUAL( -6, d_infty(t,x, 0) );
//	std::cout << "BregmanDistance d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}
