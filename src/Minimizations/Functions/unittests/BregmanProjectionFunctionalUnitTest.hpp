/*
 * BregmanProjectionFunctionalUnitTest.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef BREGMANPROJECTIONFUNCTIONALUNITTEST_HPP_
#define BREGMANPROJECTIONFUNCTIONALUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class BregmanProjectionFunctionalUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( BregmanProjectionFunctionalUnitTest );
  CPPUNIT_TEST( oneNorm );
  CPPUNIT_TEST( oneoneNorm );
  CPPUNIT_TEST( twoNorm );
  CPPUNIT_TEST( inftyNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneNorm();
  void oneoneNorm();
  void twoNorm();
  void inftyNorm();

private:
	//!> "global" tolerance threshold for minimization
	static const double tolerance;
};


#endif /* BREGMANPROJECTIONFUNCTIONALUNITTEST_HPP_ */
