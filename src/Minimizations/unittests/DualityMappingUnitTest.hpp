/*
 * DualityMappingUnitTest.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef DUALITYMAPPINGUNITTEST_HPP_
#define DUALITYMAPPINGUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class DualityMappingUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( DualityMappingUnitTest );
  CPPUNIT_TEST( oneNorm );
  CPPUNIT_TEST( twoNorm );
  CPPUNIT_TEST( inftyNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneNorm();
  void twoNorm();
  void inftyNorm();
};


#endif /* DUALITYMAPPINGUNITTEST_HPP_ */
