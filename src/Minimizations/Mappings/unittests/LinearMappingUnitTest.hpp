/*
 * LinearMappingUnitTest.hpp
 *
 *  Created on: Jul 26, 2018
 *      Author: heber
 */

#ifndef LINEARMAPPINGMAPPINGUNITTEST_HPP_
#define LINEARMAPPINGMAPPINGUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class LinearMappingUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( LinearMappingUnitTest );
  CPPUNIT_TEST( operatorTest );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void operatorTest();

private:
  //!> numeric tolerance for all tests
  static double tolerance;
};


#endif /* LINEARMAPPINGMAPPINGUNITTEST_HPP_ */
