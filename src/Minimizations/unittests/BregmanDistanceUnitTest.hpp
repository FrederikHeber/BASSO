/*
 * BregmanDistanceUnitTest.hpp
 *
 *  Created on: Apr 25, 2014
 *      Author: heber
 */

#ifndef BREGMANDISTANCEUNITTEST_HPP_
#define BREGMANDISTANCEUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class BregmanDistanceUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( BregmanDistanceUnitTest );
  CPPUNIT_TEST( oneoneNorm );
  CPPUNIT_TEST( onefiveNorm );
  CPPUNIT_TEST( twoNorm );
  CPPUNIT_TEST( threeNorm );
  CPPUNIT_TEST( sixNorm );
  CPPUNIT_TEST( inftyNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneoneNorm();
  void onefiveNorm();
  void twoNorm();
  void threeNorm();
  void sixNorm();
  void inftyNorm();
};


#endif /* BREGMANDISTANCEUNITTEST_HPP_ */
