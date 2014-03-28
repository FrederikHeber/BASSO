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


#endif /* LPNORMUNITTEST_HPP_ */
