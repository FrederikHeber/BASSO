/*
 * VectorProjectionUnitTest.hpp
 *
 *  Created on: Sep 10, 2014
 *      Author: heber
 */

#ifndef VECTORPROJECTIONUNITTEST_HPP_
#define VECTORPROJECTIONUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class VectorProjectionUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( VectorProjectionUnitTest );
  CPPUNIT_TEST( oneoneNorm );
  CPPUNIT_TEST( onefiveNorm );
  CPPUNIT_TEST( twoNorm );
  CPPUNIT_TEST( sixNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneoneNorm();
  void onefiveNorm();
  void twoNorm();
  void sixNorm();

private:

	//!> "global" tolerance threshold for minimization
	static const double tolerance;
};


#endif /* VECTORPROJECTIONUNITTEST_HPP_ */
