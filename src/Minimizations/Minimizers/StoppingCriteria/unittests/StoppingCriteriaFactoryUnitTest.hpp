/*
 * StoppingCriteriaFactoryUnitTest.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef STOPPINGCRITERIAFACTORYUNITTEST_HPP_
#define STOPPINGCRITERIAFACTORYUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>


class StoppingCriteriaFactoryUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( StoppingCriteriaFactoryUnitTest );
  CPPUNIT_TEST( singleinstanceTest );
  CPPUNIT_TEST( combinationTest );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void singleinstanceTest();
  void combinationTest();

private:
};


#endif /* STOPPINGCRITERIAFACTORYUNITTEST_HPP_ */
