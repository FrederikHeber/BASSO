/*
 * LpDualityMappingUnitTest.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef LPDUALITYMAPPINGUNITTEST_HPP_
#define LPDUALITYMAPPINGUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class LpDualityMappingUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( LpDualityMappingUnitTest );
  CPPUNIT_TEST( twoNorm );
  CPPUNIT_TEST( fourNorm );
  CPPUNIT_TEST( elevenNorm );
  CPPUNIT_TEST( otherNorm );
  CPPUNIT_TEST( setTolerance );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void twoNorm();
  void fourNorm();
  void elevenNorm();
  void otherNorm();
  void setTolerance();
};


#endif /* LPDUALITYMAPPINGUNITTEST_HPP_ */
