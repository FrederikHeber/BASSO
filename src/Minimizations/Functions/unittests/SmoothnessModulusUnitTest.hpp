/*
 * SmoothnessModulusUnitTest.hpp
 *
 *  Created on: Apr 03, 2014
 *      Author: heber
 */

#ifndef SMOOTHNESSMODULUSUNITTEST_HPP_
#define SMOOTHNESSMODULUSUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class SmoothnessModulusUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SmoothnessModulusUnitTest );
  CPPUNIT_TEST( throwTest );
  CPPUNIT_TEST( oneandhalftest );
  CPPUNIT_TEST( twotest );
  CPPUNIT_TEST( threetest );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void throwTest();
  void oneandhalftest();
  void twotest();
  void threetest();
};


#endif /* SMOOTHNESSMODULUSUNITTEST_HPP_ */
