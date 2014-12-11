/*
 * NormFactoryUnitTest.hpp
 *
 *  Created on: Dec 11, 2014
 *      Author: heber
 */

#ifndef NORMFACTORYUNITTEST_HPP_
#define NORMFACTORYUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

class NormFactoryUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( NormFactoryUnitTest );
  CPPUNIT_TEST( throwTest );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void throwTest();
};


#endif /* NORMFACTORYUNITTEST_HPP_ */
