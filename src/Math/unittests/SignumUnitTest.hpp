/*
 * SignumUnitTest.hpp
 *
 *  Created on: Jul 14, 2014
 *      Author: heber
 */

#ifndef SIGNUMUNITTEST_HPP_
#define SIGNUMUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class SignumUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SignumUnitTest );
  CPPUNIT_TEST( twoTest );
  CPPUNIT_TEST( fiveTest );
  CPPUNIT_TEST( tenTest );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void twoTest();
  void fiveTest();
  void tenTest();
};


#endif /* SIGNUMUNITTEST_HPP_ */
