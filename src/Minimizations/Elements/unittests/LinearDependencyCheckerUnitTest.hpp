/*
 * LinearDependencyCheckerUnitTest.hpp
 *
 *  Created on: Feb 10, 2015
 *      Author: heber
 */

#ifndef LINEARDEPENDENCYCHECKERUNITTEST_HPP_
#define LINEARDEPENDENCYCHECKERUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

#include "Minimizations/Elements/LinearDependencyChecker.hpp"

class LinearDependencyCheckerUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( LinearDependencyCheckerUnitTest );
  CPPUNIT_TEST( noVector );
  CPPUNIT_TEST( singleVector );
  CPPUNIT_TEST( twoVectors );
  CPPUNIT_TEST( threeVectors );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void noVector();
  void singleVector();
  void twoVectors();
  void threeVectors();

private:
  LinearDependencyChecker lindepcheck;
};


#endif /* LINEARDEPENDENCYCHECKERUNITTEST_HPP_ */
