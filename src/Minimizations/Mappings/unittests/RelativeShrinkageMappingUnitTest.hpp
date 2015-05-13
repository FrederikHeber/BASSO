/*
 * RelativeShrinkageMappingUnitTest.hpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#ifndef SOFTTHRESHOLDINGMAPPINGUNITTEST_HPP_
#define SOFTTHRESHOLDINGMAPPINGUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class RelativeShrinkageMappingUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( RelativeShrinkageMappingUnitTest );
  CPPUNIT_TEST( oneNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneNorm();

private:
  //!> numeric tolerance for all tests
  static double tolerance;
};


#endif /* SOFTTHRESHOLDINGMAPPINGUNITTEST_HPP_ */
