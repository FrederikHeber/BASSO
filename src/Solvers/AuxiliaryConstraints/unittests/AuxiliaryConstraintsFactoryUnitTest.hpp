/*
 * AuxiliaryConstraintsFactoryUnitTest.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef AUXILIARYCONSTRAINTSFACTORYUNITTEST_HPP_
#define AUXILIARYCONSTRAINTSFACTORYUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>


class AuxiliaryConstraintsFactoryUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( AuxiliaryConstraintsFactoryUnitTest );
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


#endif /* AUXILIARYCONSTRAINTSFACTORYUNITTEST_HPP_ */
