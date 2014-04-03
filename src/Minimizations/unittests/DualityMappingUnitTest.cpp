/*
 * DualityMappingUnitTest.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "DualityMappingUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/DualityMapping.hpp"
#include "Minimizations/MinimizationExceptions.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DualityMappingUnitTest );


void DualityMappingUnitTest::setUp()
{
}


void DualityMappingUnitTest::tearDown()
{
}

void DualityMappingUnitTest::throwTest()
{
	// we check that assertion is thrown for invalid p value
	std::cout << "The following assertion is intended and does not indicate a failure of the test." << std::endl;
	CPPUNIT_ASSERT_THROW(
			DualityMapping J_illegal(-0.5),
			MinimizationIllegalValue_exception );
}

void DualityMappingUnitTest::oneNorm()
{
	Eigen::VectorXd v(2);
	DualityMapping J_1(1);
	v << 4,3;
	{
		Eigen::VectorXd compare(2);
		compare << 1, 1;
//		std::cout << "DualityMapping J_1 with weight 1 of v is ("
//				<< J_1(v,1).transpose() << ")" << std::endl;
		CPPUNIT_ASSERT( J_1(v,1).isApprox(compare) );
	}
	{
		Eigen::VectorXd compare(2);
		compare << 7, 7;
//		std::cout << "DualityMapping J_1 with weight 2 of v is ("
//				<< J_1(v,2).transpose() << ")" << std::endl;
		CPPUNIT_ASSERT( J_1(v,2).isApprox(compare) );
	}
	{
		Eigen::VectorXd compare(2);
		compare << 49, 49;
//		std::cout << "DualityMapping J_1 with weight 3 of v is ("
//				<< J_1(v,3).transpose() << ")" << std::endl;
		CPPUNIT_ASSERT( J_1(v,3).isApprox(compare) );
	}
}

void DualityMappingUnitTest::twoNorm()
{
	Eigen::VectorXd v(2);
	DualityMapping J_2(2);
	v << 4,3;
	{
		Eigen::VectorXd compare(2);
		compare << 4, 3;
//		std::cout << "DualityMapping J_2 with weight 2 of v is ("
//				<< J_2(v,2).transpose() << ")" << std::endl;
		CPPUNIT_ASSERT( J_2(v,2).isApprox(compare) );
	}
	{
		Eigen::VectorXd compare(2);
		compare << 0.8, 0.6;
//		std::cout << "DualityMapping J_2 with weight 1 of v is ("
//				<< J_2(v,1).transpose() << ")" << std::endl;
		CPPUNIT_ASSERT( J_2(v,1).isApprox(compare) );
	}
	{
		Eigen::VectorXd compare(2);
		compare << 100, 75;
//		std::cout << "DualityMapping J_2 with weight 4 of v is ("
//				<< J_2(v,4).transpose() << ")" << std::endl;
		CPPUNIT_ASSERT( J_2(v,4).isApprox(compare) );
	}
}

void DualityMappingUnitTest::inftyNorm()
{
	Eigen::VectorXd v(2);
	DualityMapping J_infty(0);
	v << 4,3;
	{
		Eigen::VectorXd compare(2);
		compare << 1, 0;
//		std::cout << "DualityMapping J_infty with weight 1 of v is ("
//				<< J_infty(v,1).transpose() << ")" << std::endl;
		CPPUNIT_ASSERT( J_infty(v,1).isApprox(compare) );
	}
	{
		Eigen::VectorXd compare(2);
		compare << 4, 0;
//		std::cout << "DualityMapping J_infty with weight 2 of v is ("
//				<< J_infty(v,2).transpose() << ")" << std::endl;
		CPPUNIT_ASSERT( J_infty(v,2).isApprox(compare) );
	}
	{
		Eigen::VectorXd compare(2);
		compare << 16, 0;
//		std::cout << "DualityMapping J_infty with weight 3 of v is ("
//				<< J_infty(v,3).transpose() << ")" << std::endl;
		CPPUNIT_ASSERT( J_infty(v,3).isApprox(compare) );
	}
}

void DualityMappingUnitTest::setTolerance()
{
	DualityMapping J_infty(0);
	CPPUNIT_ASSERT_EQUAL(BASSOTOLERANCE, J_infty.tolerance);
	const double value = 1e-1;
	J_infty.setTolerance(value);
	CPPUNIT_ASSERT_EQUAL(value, J_infty.tolerance);
}
