/*
 * point_tUnitTest.hpp
 *
 *  Created on: Jul 07, 2015
 *      Author: heber
 */

#ifndef POINT_TUNITTEST_HPP_
#define POINT_TUNITTEST_HPP_

#include "BassoConfig.h"

#include <cppunit/extensions/HelperMacros.h>


class point_tUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( point_tUnitTest );
  CPPUNIT_TEST( getSquaredDistanceTest );
  CPPUNIT_TEST( maxNormTest );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void getSquaredDistanceTest();
  void maxNormTest();
};


#endif /* POINT_TUNITTEST_HPP_ */
