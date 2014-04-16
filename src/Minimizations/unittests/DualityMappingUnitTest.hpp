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
  CPPUNIT_TEST( throwTest );
  CPPUNIT_TEST( oneNorm );
  CPPUNIT_TEST( twoNorm );
  CPPUNIT_TEST( fourNorm );
  CPPUNIT_TEST( elevenNorm );
  CPPUNIT_TEST( inftyNorm );
  CPPUNIT_TEST( setTolerance );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void throwTest();
  void oneNorm();
  void twoNorm();
  void fourNorm();
  void elevenNorm();
  void inftyNorm();
  void setTolerance();
};


#endif /* DUALITYMAPPINGUNITTEST_HPP_ */
