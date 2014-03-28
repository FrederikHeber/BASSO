/*
 * BregmanFunctionalUnitTest.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BregmanFunctionalUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/BregmanFunctional.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( BregmanFunctionalUnitTest );


void BregmanFunctionalUnitTest::setUp()
{
}


void BregmanFunctionalUnitTest::tearDown()
{
}

void BregmanFunctionalUnitTest::oneNorm()
{
	BregmanFunctional d_1(1);
	Eigen::VectorXd t(2);
	t << 4,3;
	Eigen::VectorXd x(2);
	x << 4,3;
	Eigen::MatrixXd U(2,2);
	U << 1,0,0,1;
	Eigen::VectorXd alpha(2);
	alpha << 1,0;
	const unsigned int q = 2; 		// power of weight of duality mapping
	Eigen::VectorXd compare(2);
	compare << 1, 0;
	CPPUNIT_ASSERT_EQUAL( 4., d_1(t,x,U,alpha,q) );
	CPPUNIT_ASSERT( d_1.gradient(t,x,U,alpha,q).isApprox(compare) );
//	std::cout << "BregmanFunctional d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}

void BregmanFunctionalUnitTest::twoNorm()
{
	BregmanFunctional d_2(2);
	Eigen::VectorXd t(2);
	t << 4,3;
	Eigen::VectorXd x(2);
	x << 4,3;
	Eigen::MatrixXd U(2,2);
	U << 1,0,0,1;
	Eigen::VectorXd alpha(2);
	alpha << 1,0;
	const unsigned int q = 2; 		// power of weight of duality mapping
	Eigen::VectorXd compare(2);
	compare << 1, 0;
	CPPUNIT_ASSERT_EQUAL( 4., d_2(t,x,U,alpha,q) );
	CPPUNIT_ASSERT( d_2.gradient(t,x,U,alpha,q).isApprox(compare) );
//	std::cout << "BregmanFunctional d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}

void BregmanFunctionalUnitTest::inftyNorm()
{
	BregmanFunctional d_infty(LpNorm::Infinity);
	Eigen::VectorXd t(2);
	t << 4,3;
	Eigen::VectorXd x(2);
	x << 4,3;
	Eigen::MatrixXd U(2,2);
	U << 1,0,0,1;
	Eigen::VectorXd alpha(2);
	alpha << 1,0;
	const unsigned int q = 2; 		// power of weight of duality mapping
	Eigen::VectorXd compare(2);
	compare << 1, 0;
	CPPUNIT_ASSERT_EQUAL( 4., d_infty(t,x,U,alpha,q) );
	CPPUNIT_ASSERT( d_infty.gradient(t,x,U,alpha,q).isApprox(compare) );
//	std::cout << "BregmanFunctional d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}
