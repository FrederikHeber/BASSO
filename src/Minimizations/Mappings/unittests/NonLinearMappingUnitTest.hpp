/*
 * NonLinearMappingUnitTest.hpp
 *
 *  Created on: Jul 26, 2018
 *      Author: heber
 */

#ifndef NONLINEARMAPPINGMAPPINGUNITTEST_HPP_
#define NONLINEARMAPPINGMAPPINGUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class NonLinearMappingUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( NonLinearMappingUnitTest );
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


#endif /* NONLINEARMAPPINGMAPPINGUNITTEST_HPP_ */
