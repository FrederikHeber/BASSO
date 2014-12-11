/*
 * LpNormUnitTest.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef LPNORMUNITTEST_HPP_
#define LPNORMUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class LpNormUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( LpNormUnitTest );
  CPPUNIT_TEST( twoNorm );
  CPPUNIT_TEST( threeNorm );
  CPPUNIT_TEST( fourNorm );
  CPPUNIT_TEST( elevenNorm );
  CPPUNIT_TEST( otherNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void twoNorm();
  void threeNorm();
  void fourNorm();
  void elevenNorm();
  void otherNorm();
};


#endif /* LPNORMUNITTEST_HPP_ */
