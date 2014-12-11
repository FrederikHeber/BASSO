/*
 * LInfinityDualityMappingUnitTest.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef LINFINITYDUALITYMAPPINGUNITTEST_HPP_
#define LINFINITYDUALITYMAPPINGUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class LInfinityDualityMappingUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( LInfinityDualityMappingUnitTest );
  CPPUNIT_TEST( inftyNorm );
  CPPUNIT_TEST( setTolerance );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void inftyNorm();
  void setTolerance();

private:
  //!> numeric tolerance for all tests
  static double tolerance;
};


#endif /* LINFINITYDUALITYMAPPINGUNITTEST_HPP_ */
