/*
 * L1DualityMappingUnitTest.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef L1DUALITYMAPPINGUNITTEST_HPP_
#define L1DUALITYMAPPINGUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class L1DualityMappingUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( L1DualityMappingUnitTest );
  CPPUNIT_TEST( oneNorm );
  CPPUNIT_TEST( setTolerance );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneNorm();
  void setTolerance();
};


#endif /* L1DUALITYMAPPINGUNITTEST_HPP_ */
