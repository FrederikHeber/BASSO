/*
 * IllegalDualityMappingUnitTest.cpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "IllegalDualityMappingUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/DualityMappings/IllegalDualityMapping.hpp"
#include "Minimizations/MinimizationExceptions.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( IllegalDualityMappingUnitTest );


void IllegalDualityMappingUnitTest::setUp()
{
}


void IllegalDualityMappingUnitTest::tearDown()
{
}

void IllegalDualityMappingUnitTest::IllegalCall()
{
	IllegalDualityMapping J_illegal;
	Eigen::VectorXd x(10);
	x << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929;
	CPPUNIT_ASSERT_THROW( J_illegal(x, .9), MinimizationIllegalValue_exception);
}
