/*
 * SoftThresholdingOperatorUnitTest.hpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#ifndef SOFTTHRESHOLDINGOPERATORUNITTEST_HPP_
#define SOFTTHRESHOLDINGOPERATORUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class SoftThresholdingOperatorUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SoftThresholdingOperatorUnitTest );
  CPPUNIT_TEST( oneNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneNorm();
};


#endif /* SOFTTHRESHOLDINGOPERATORUNITTEST_HPP_ */
