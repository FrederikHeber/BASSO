/*
 * StoppingCriteriaFactoryUnitTest.cpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#include <boost/assign.hpp>
#include <iostream>

#include "Minimizations/Minimizers/StoppingCriteria/unittests/StoppingCriteriaFactoryUnitTest.hpp"

#include "Log/Logging.hpp"

#include "Minimizations/Minimizers/StoppingCriteria/CheckDivergentResiduum.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/CheckIterationCount.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/CheckRelativeResiduum.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/CheckResiduum.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/CheckWalltime.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_AND.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_NOT.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriterion_OR.hpp"
#include "Minimizations/Minimizers/StoppingCriteria/StoppingCriteriaFactory.hpp"

using namespace boost::assign;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( StoppingCriteriaFactoryUnitTest );

void StoppingCriteriaFactoryUnitTest::setUp()
{
	// StoppingCriteriaFactory uses logging
	boost::log::core::get()->set_filter
			(
					boost::log::trivial::severity >= boost::log::trivial::info
			);
	startLogging();
}


void StoppingCriteriaFactoryUnitTest::tearDown()
{
}

void StoppingCriteriaFactoryUnitTest::singleinstanceTest()
{
	StoppingCriteriaFactory factory;
	StoppingArguments args;

	{
		std::string criteria_line = "DivergentResiduum";
		StoppingCriterion::ptr_t criterion =
				factory.create(criteria_line, args);
		CPPUNIT_ASSERT( dynamic_cast<CheckDivergentResiduum *>(criterion.get()) != NULL);
	}
	{
		std::string criteria_line = "MaxIterationCount";
		StoppingCriterion::ptr_t criterion =
				factory.create(criteria_line, args);
		CPPUNIT_ASSERT( dynamic_cast<CheckIterationCount *>(criterion.get()) != NULL);
	}
	{
		std::string criteria_line = "RelativeResiduum";
		StoppingCriterion::ptr_t criterion =
				factory.create(criteria_line, args);
		CPPUNIT_ASSERT( dynamic_cast<CheckRelativeResiduum *>(criterion.get()) != NULL);
	}
	{
		std::string criteria_line = "Residuum";
		StoppingCriterion::ptr_t criterion =
				factory.create(criteria_line, args);
		CPPUNIT_ASSERT( dynamic_cast<CheckResiduum *>(criterion.get()) != NULL);
	}
	{
		std::string criteria_line = "MaxWalltime";
		StoppingCriterion::ptr_t criterion =
				factory.create(criteria_line, args);
		CPPUNIT_ASSERT( dynamic_cast<CheckWalltime *>(criterion.get()) != NULL);
	}
}

void StoppingCriteriaFactoryUnitTest::combinationTest()
{
	StoppingCriteriaFactory factory;
	StoppingArguments args;

	// AND
	{
		std::string criteria_line = "MaxIterationCount && Residuum";
		StoppingCriterion::ptr_t criterion =
				factory.create(criteria_line, args);
		CPPUNIT_ASSERT( dynamic_cast<StoppingCriterion_AND *>(criterion.get()) != NULL);
		StoppingCriterion_AND &combined =
				dynamic_cast<StoppingCriterion_AND &>(*criterion.get());
		CPPUNIT_ASSERT( dynamic_cast<CheckIterationCount *>(combined.left.get()) != NULL);
		CPPUNIT_ASSERT( dynamic_cast<CheckResiduum *>(combined.right.get()) != NULL);
	}

	// NOT
	{
		std::string criteria_line = "! MaxIterationCount";
		StoppingCriterion::ptr_t criterion =
				factory.create(criteria_line, args);
		CPPUNIT_ASSERT( dynamic_cast<StoppingCriterion_NOT *>(criterion.get()) != NULL);
		StoppingCriterion_NOT &combined =
				dynamic_cast<StoppingCriterion_NOT &>(*criterion.get());
		CPPUNIT_ASSERT( dynamic_cast<CheckIterationCount *>(combined.impl.get()) != NULL);
	}

	// OR
	{
		std::string criteria_line = "MaxIterationCount || Residuum";
		StoppingCriterion::ptr_t criterion =
				factory.create(criteria_line, args);
		CPPUNIT_ASSERT( dynamic_cast<StoppingCriterion_OR *>(criterion.get()) != NULL);
		StoppingCriterion_OR &combined =
				dynamic_cast<StoppingCriterion_OR &>(*criterion.get());
		CPPUNIT_ASSERT( dynamic_cast<CheckIterationCount *>(combined.left.get()) != NULL);
		CPPUNIT_ASSERT( dynamic_cast<CheckResiduum *>(combined.right.get()) != NULL);
	}

}
