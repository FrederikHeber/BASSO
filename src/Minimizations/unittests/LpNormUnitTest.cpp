/*
 * LpNormUnitTest.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "LpNormUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/LpNorm.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( LpNormUnitTest );


void LpNormUnitTest::setUp()
{
}


void LpNormUnitTest::tearDown()
{
}

void LpNormUnitTest::oneNorm()
{
	LpNorm Norm1(1);
	{
		Eigen::VectorXd x(2);
		x << 2,1;
		const double expected = x.lpNorm<1>();
		CPPUNIT_ASSERT_EQUAL( expected, Norm1(x));
	}
	{
		Eigen::VectorXd x(2);
		x << 1,0;
		const double expected = x.lpNorm<1>();
		CPPUNIT_ASSERT_EQUAL( expected, Norm1(x));
	}
	{
		Eigen::VectorXd x(2);
		x << 1,1;
		const double expected = x.lpNorm<1>();
		CPPUNIT_ASSERT_EQUAL( expected, Norm1(x));
	}
	{
		Eigen::VectorXd x(2);
		x << -2,1;
		const double expected = x.lpNorm<1>();
		CPPUNIT_ASSERT_EQUAL( expected, Norm1(x));
	}
	{
		Eigen::VectorXd x(2);
		x << 2,-1;
		const double expected = x.lpNorm<1>();
		CPPUNIT_ASSERT_EQUAL( expected, Norm1(x));
	}
}

void LpNormUnitTest::twoNorm()
{
	LpNorm Norm2(2);
	{
		Eigen::VectorXd x(2);
		x << 2,1;
		const double expected = x.lpNorm<2>();
		CPPUNIT_ASSERT_EQUAL( expected, Norm2(x));
	}
	{
		Eigen::VectorXd x(2);
		x << 1,0;
		const double expected = x.lpNorm<2>();
		CPPUNIT_ASSERT_EQUAL( expected, Norm2(x));
	}
	{
		Eigen::VectorXd x(2);
		x << 1,1;
		const double expected = x.lpNorm<2>();
		CPPUNIT_ASSERT_EQUAL( expected, Norm2(x));
	}
	{
		Eigen::VectorXd x(2);
		x << -2,1;
		const double expected = x.lpNorm<2>();
		CPPUNIT_ASSERT_EQUAL( expected, Norm2(x));
	}
	{
		Eigen::VectorXd x(2);
		x << 2,-1;
		const double expected = x.lpNorm<2>();
		CPPUNIT_ASSERT_EQUAL( expected, Norm2(x));
	}
}

void LpNormUnitTest::inftyNorm()
{
	LpNorm NormInfty(0);
	{
		Eigen::VectorXd x(2);
		x << 2,1;
		const double expected = x.lpNorm<Eigen::Infinity>();
		CPPUNIT_ASSERT_EQUAL( expected, NormInfty(x));
	}
	{
		Eigen::VectorXd x(2);
		x << 1,0;
		const double expected = x.lpNorm<Eigen::Infinity>();
		CPPUNIT_ASSERT_EQUAL( expected, NormInfty(x));
	}
	{
		Eigen::VectorXd x(2);
		x << 1,1;
		const double expected = x.lpNorm<Eigen::Infinity>();
		CPPUNIT_ASSERT_EQUAL( expected, NormInfty(x));
	}
	{
		Eigen::VectorXd x(2);
		x << -2,1;
		const double expected = x.lpNorm<Eigen::Infinity>();
		CPPUNIT_ASSERT_EQUAL( expected, NormInfty(x));
	}
	{
		Eigen::VectorXd x(2);
		x << 2,-1;
		const double expected = x.lpNorm<Eigen::Infinity>();
		CPPUNIT_ASSERT_EQUAL( expected, NormInfty(x));
	}
}