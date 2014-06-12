/*
 * DatabaseUnitTest.hpp
 *
 *  Created on: Jun 12, 2014
 *      Author: heber
 */

#ifndef DATABASEUNITTEST_HPP_
#define DATABASEUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class DatabaseUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( DatabaseUnitTest );
  CPPUNIT_TEST( simpleTest );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void simpleTest();
};


#endif /* DATABASEUNITTEST_HPP_ */
