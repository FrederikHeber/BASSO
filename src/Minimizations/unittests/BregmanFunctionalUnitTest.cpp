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
	BregmanFunctional<1> d_1;
	Eigen::VectorXd t(2);
	t << 4,3;
	Eigen::VectorXd x(2);
	x << 4,3;
	Eigen::MatrixXd U(2,2);
	U << 1,0,0,1;
	Eigen::VectorXd alpha(2);
	alpha << 1,0;
	const unsigned int q = 2; 		// power of weight of duality mapping
	std::pair<double, Eigen::VectorXd> vals = d_1(t,x,U,alpha,q);
	Eigen::VectorXd compare(2);
	compare << 1, 0;
	CPPUNIT_ASSERT_EQUAL( 4., vals.first );
	CPPUNIT_ASSERT( vals.second.isApprox(compare) );
//	std::cout << "BregmanFunctional d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}

void BregmanFunctionalUnitTest::twoNorm()
{
	BregmanFunctional<2> d_2;
	Eigen::VectorXd t(2);
	t << 4,3;
	Eigen::VectorXd x(2);
	x << 4,3;
	Eigen::MatrixXd U(2,2);
	U << 1,0,0,1;
	Eigen::VectorXd alpha(2);
	alpha << 1,0;
	const unsigned int q = 2; 		// power of weight of duality mapping
	std::pair<double, Eigen::VectorXd> vals = d_2(t,x,U,alpha,q);
	Eigen::VectorXd compare(2);
	compare << 1, 0;
	CPPUNIT_ASSERT_EQUAL( 4., vals.first );
	CPPUNIT_ASSERT( vals.second.isApprox(compare) );
//	std::cout << "BregmanFunctional d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}

void BregmanFunctionalUnitTest::inftyNorm()
{
	BregmanFunctional<Eigen::Infinity> d_infty;
	Eigen::VectorXd t(2);
	t << 4,3;
	Eigen::VectorXd x(2);
	x << 4,3;
	Eigen::MatrixXd U(2,2);
	U << 1,0,0,1;
	Eigen::VectorXd alpha(2);
	alpha << 1,0;
	const unsigned int q = 2; 		// power of weight of duality mapping
	std::pair<double, Eigen::VectorXd> vals = d_infty(t,x,U,alpha,q);
	Eigen::VectorXd compare(2);
	compare << 1, 0;
	CPPUNIT_ASSERT_EQUAL( 4., vals.first );
	CPPUNIT_ASSERT( vals.second.isApprox(compare) );
//	std::cout << "BregmanFunctional d_2 of v is "
//			<< vals.first << "," << vals.second.transpose() << "" << std::endl;
}
