/*
 * BregmanFunctionalUnitTest.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef BREGMANFUNCTIONALUNITTEST_HPP_
#define BREGMANFUNCTIONALUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class BregmanFunctionalUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( BregmanFunctionalUnitTest );
  CPPUNIT_TEST( oneoneNorm );
  CPPUNIT_TEST( twoNorm );
  CPPUNIT_TEST( inftyNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneoneNorm();
  void twoNorm();
  void inftyNorm();
};


#endif /* BREGMANFUNCTIONALUNITTEST_HPP_ */
