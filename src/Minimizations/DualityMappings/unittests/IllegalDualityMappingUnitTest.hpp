/*
 * IllegalDualityMappingUnitTest.hpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#ifndef ILLEGALDUALITYMAPPINGUNITTEST_HPP_
#define ILLEGALDUALITYMAPPINGUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class IllegalDualityMappingUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( IllegalDualityMappingUnitTest );
  CPPUNIT_TEST( IllegalCall );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void IllegalCall();
};


#endif /* ILLEGALDUALITYMAPPINGUNITTEST_HPP_ */
