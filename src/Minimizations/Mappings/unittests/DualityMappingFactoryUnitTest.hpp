/*
 * DualityMappingFactoryUnitTest.hpp
 *
 *  Created on: Jul 01, 2015
 *      Author: heber
 */

#ifndef DUALITYMAPPINGFACTORYUNITTEST_HPP_
#define DUALITYMAPPINGFACTORYUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class DualityMappingFactoryUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( DualityMappingFactoryUnitTest );
  CPPUNIT_TEST( equivalentTokensTest );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void equivalentTokensTest();
};


#endif /* DUALITYMAPPINGFACTORYUNITTEST_HPP_ */
