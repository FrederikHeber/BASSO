/*
 * L1NormUnitTest.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef L1NORMUNITTEST_HPP_
#define L1NORMUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class L1NormUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( L1NormUnitTest );
  CPPUNIT_TEST( oneNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneNorm();
};


#endif /* L1NORMUNITTEST_HPP_ */
