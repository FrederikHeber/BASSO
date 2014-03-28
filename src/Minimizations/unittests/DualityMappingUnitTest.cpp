/*
 * DualityMappingUnitTest.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "DualityMappingUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/DualityMapping.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DualityMappingUnitTest );


void DualityMappingUnitTest::setUp()
{
}


void DualityMappingUnitTest::tearDown()
{
}


void DualityMappingUnitTest::oneNorm()
{
	Eigen::VectorXd v(2);
	v << 4,3;
	{
		DualityMapping<1> J_1(1);
		Eigen::VectorXd compare(2);
		compare << 1, 1;
		CPPUNIT_ASSERT( J_1(v).isApprox(compare) );
//		std::cout << "DualityMapping J_1 with weight 1 of v is ("
//				<< J_1(v).transpose() << ")" << std::endl;
	}
	{
		DualityMapping<1> J_1(2);
		Eigen::VectorXd compare(2);
		compare << 7, 7;
		CPPUNIT_ASSERT( J_1(v).isApprox(compare) );
//		std::cout << "DualityMapping J_1 with weight 2 of v is ("
//				<< J_1(v).transpose() << ")" << std::endl;
	}
	{
		DualityMapping<1> J_1(3);
		Eigen::VectorXd compare(2);
		compare << 49, 49;
		CPPUNIT_ASSERT( J_1(v).isApprox(compare) );
//		std::cout << "DualityMapping J_1 with weight 3 of v is ("
//				<< J_1(v).transpose() << ")" << std::endl;
	}
}

void DualityMappingUnitTest::twoNorm()
{
	Eigen::VectorXd v(2);
	v << 4,3;
	{
		DualityMapping<2> J_2(2);
		Eigen::VectorXd compare(2);
		compare << 4, 3;
		CPPUNIT_ASSERT( J_2(v).isApprox(compare) );
//		std::cout << "DualityMapping J_2 with weight 2 of v is ("
//				<< J_2(v).transpose() << ")" << std::endl;
	}
	{
		DualityMapping<2> J_2(1);
		Eigen::VectorXd compare(2);
		compare << 0.8, 0.6;
		CPPUNIT_ASSERT( J_2(v).isApprox(compare) );
//		std::cout << "DualityMapping J_2 with weight 1 of v is ("
//				<< J_2(v).transpose() << ")" << std::endl;
	}
	{
		DualityMapping<2> J_2(4);
		Eigen::VectorXd compare(2);
		compare << 100, 75;
		CPPUNIT_ASSERT( J_2(v).isApprox(compare) );
//		std::cout << "DualityMapping J_2 with weight 4 of v is ("
//				<< J_2(v).transpose() << ")" << std::endl;
	}
}

void DualityMappingUnitTest::inftyNorm()
{
	Eigen::VectorXd v(2);
	v << 4,3;
	{
		DualityMapping<Eigen::Infinity> J_infty(1);
		Eigen::VectorXd compare(2);
		compare << 1, 0;
		CPPUNIT_ASSERT( J_infty(v).isApprox(compare) );
//		std::cout << "DualityMapping J_infty with weight 1 of v is ("
//				<< J_infty(v).transpose() << ")" << std::endl;
	}
	{
		DualityMapping<Eigen::Infinity> J_infty(2);
		Eigen::VectorXd compare(2);
		compare << 4, 0;
		CPPUNIT_ASSERT( J_infty(v).isApprox(compare) );
//		std::cout << "DualityMapping J_infty with weight 2 of v is ("
//				<< J_infty(v).transpose() << ")" << std::endl;
	}
	{
		DualityMapping<Eigen::Infinity> J_infty(3);
		Eigen::VectorXd compare(2);
		compare << 16, 0;
		CPPUNIT_ASSERT( J_infty(v).isApprox(compare) );
//		std::cout << "DualityMapping J_infty with weight 3 of v is ("
//				<< J_infty(v).transpose() << ")" << std::endl;
	}
}

