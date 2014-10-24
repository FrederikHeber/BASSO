/*
 * LInfinityNormUnitTest.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef LINFINITYNORMUNITTEST_HPP_
#define LINFINITYNORMUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class LInfinityNormUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( LInfinityNormUnitTest );
  CPPUNIT_TEST( inftyNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void inftyNorm();
};


#endif /* LINFINITYNORMUNITTEST_HPP_ */
