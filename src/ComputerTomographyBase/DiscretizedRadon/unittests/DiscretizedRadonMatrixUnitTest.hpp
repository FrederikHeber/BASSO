/*
 * DiscretizedRadonMatrixUnitTest.hpp
 *
 *  Created on: Jul 07, 2015
 *      Author: heber
 */

#ifndef DISCRETIZEDRADONMATRIXUNITTEST_HPP_
#define DISCRETIZEDRADONMATRIXUNITTEST_HPP_

#include "BassoConfig.h"

#include <cppunit/extensions/HelperMacros.h>


class DiscretizedRadonMatrixUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( DiscretizedRadonMatrixUnitTest );
  CPPUNIT_TEST( adjointTest );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void adjointTest();
};


#endif /* DISCRETIZEDRADONMATRIXUNITTEST_HPP_ */
