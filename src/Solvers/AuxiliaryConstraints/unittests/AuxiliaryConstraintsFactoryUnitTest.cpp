/*
 * AuxiliaryConstraintsFactoryUnitTest.cpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#include <boost/assign.hpp>
#include <iostream>

#include "Solvers/AuxiliaryConstraints/unittests/AuxiliaryConstraintsFactoryUnitTest.hpp"

#include "Log/Logging.hpp"

#include "Solvers/AuxiliaryConstraints/NonnegativeConstraint.hpp"
#include "Solvers/AuxiliaryConstraints/NonpositiveConstraint.hpp"
#include "Solvers/AuxiliaryConstraints/UnityConstraint.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints_AND.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraintsFactory.hpp"

using namespace boost::assign;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( AuxiliaryConstraintsFactoryUnitTest );

void AuxiliaryConstraintsFactoryUnitTest::setUp()
{
	// AuxiliaryConstraintsFactory uses logging
	boost::log::core::get()->set_filter
			(
					boost::log::trivial::severity >= boost::log::trivial::info
			);
	startLogging();
}


void AuxiliaryConstraintsFactoryUnitTest::tearDown()
{
}

void AuxiliaryConstraintsFactoryUnitTest::singleinstanceTest()
{
	AuxiliaryConstraintsFactory factory;

	{
		std::string criteria_line = "Nonnegative";
		AuxiliaryConstraints::ptr_t criterion =
				factory.create(criteria_line);
		CPPUNIT_ASSERT( dynamic_cast<NonnegativeConstraint *>(criterion.get()) != NULL);
	}
	{
		std::string criteria_line = "Nonpositive";
		AuxiliaryConstraints::ptr_t criterion =
				factory.create(criteria_line);
		CPPUNIT_ASSERT( dynamic_cast<NonpositiveConstraint *>(criterion.get()) != NULL);
	}
	{
		std::string criteria_line = "Unity";
		AuxiliaryConstraints::ptr_t criterion =
				factory.create(criteria_line);
		CPPUNIT_ASSERT( dynamic_cast<UnityConstraint *>(criterion.get()) != NULL);
	}
}

void AuxiliaryConstraintsFactoryUnitTest::combinationTest()
{
	AuxiliaryConstraintsFactory factory;

	// AND
	{
		std::string criteria_line = "Nonnegative && Unity";
		AuxiliaryConstraints::ptr_t criterion =
				factory.create(criteria_line);
		CPPUNIT_ASSERT( dynamic_cast<AuxiliaryConstraints_AND *>(criterion.get()) != NULL);
		AuxiliaryConstraints_AND &combined =
				dynamic_cast<AuxiliaryConstraints_AND &>(*criterion.get());
		CPPUNIT_ASSERT( dynamic_cast<NonnegativeConstraint *>(combined.left.get()) != NULL);
		CPPUNIT_ASSERT( dynamic_cast<UnityConstraint *>(combined.right.get()) != NULL);
	}
}
