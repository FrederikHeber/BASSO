/*
 * MetricProjectionFunctionalUnitTest.hpp
 *
 *  Created on: May 20, 2015
 *      Author: heber
 */

#ifndef METRICPROJECTIONFUNCTIONALUNITTEST_HPP_
#define METRICPROJECTIONFUNCTIONALUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class MetricProjectionFunctionalUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MetricProjectionFunctionalUnitTest );
//  CPPUNIT_TEST( oneNorm );
//  CPPUNIT_TEST( oneoneNorm );
  CPPUNIT_TEST( twoNorm );
//  CPPUNIT_TEST( inftyNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneNorm();
  void oneoneNorm();
  void twoNorm();
  void inftyNorm();

private:
	//!> "global" tolerance threshold for minimization
	static const double tolerance;
};


#endif /* METRICPROJECTIONFUNCTIONALUNITTEST_HPP_ */
