/*
 * ConjugateValueUnitTest.hpp
 *
 *  Created on: Oct 28, 2014
 *      Author: heber
 */

#ifndef CONJUGATEVALUEUNITTEST_HPP_
#define CONJUGATEVALUEUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class ConjugateValueUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ConjugateValueUnitTest );
  CPPUNIT_TEST( Test );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void Test();
};


#endif /* CONJUGATEVALUEUNITTEST_HPP_ */
