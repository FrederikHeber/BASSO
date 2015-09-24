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
  CPPUNIT_TEST( singlePixel_few_angles_few_offsetsTest );
  CPPUNIT_TEST( singlePixel_many_angles_few_offsetsTest );
  CPPUNIT_TEST( singlePixel_few_angles_many_offsetsTest );
  CPPUNIT_TEST( ninePixel_few_angles_few_offsetsTest );
  CPPUNIT_TEST( ninePixel_few_angles_many_offsetsTest );
  CPPUNIT_TEST( adjointTest );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void singlePixel_few_angles_few_offsetsTest();
  void singlePixel_many_angles_few_offsetsTest();
  void singlePixel_few_angles_many_offsetsTest();
  void ninePixel_few_angles_few_offsetsTest();
  void ninePixel_few_angles_many_offsetsTest();
  void adjointTest();
};


#endif /* DISCRETIZEDRADONMATRIXUNITTEST_HPP_ */
