/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2018 Frederik Heber. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*
 * LinearMappingUnitTest.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LinearMappingUnitTest.hpp"

#include <boost/assign.hpp>

#include <Eigen/Dense>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"

#include "../MappingFactory.hpp"
#include "Minimizations/Norms/NormExceptions.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( LinearMappingUnitTest );

using namespace boost::assign;

// static entities
double LinearMappingUnitTest::tolerance = 1e-4;


void LinearMappingUnitTest::setUp()
{
}


void LinearMappingUnitTest::tearDown()
{
}

/** We generate test vectors as follows via octave:
 *
 * x=ones(2,10)-2.*rand(2,10)
 * gval=zeros(2,10)
 * for i=1:10
 * 	gval(:,i)=LinearMapping(x, p, power, 1e-6)
 * endfor
 * gval
 *
 *
 */

void LinearMappingUnitTest::operatorTest()
{
	const double p = 2.;
	Eigen::MatrixXd matrix(2,10);
	matrix << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929,
			0.608996,0.567817,0.261505,-0.294843,-0.387682,-0.513624,-0.728372,-0.676635,-0.503763,-0.611381;

	// spaces are irrelevant here
	const double power = 2.;
	NormedSpace_ptr_t SpaceX;
	NormedSpace_ptr_t SpaceY;
	{
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		SpaceX = NormedSpaceFactory::create(
				matrix.outerSize(), "lp", args);
	}
	{
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		SpaceY = NormedSpaceFactory::create(
				matrix.innerSize(), "lp", args);
	}
	{
		Mapping_ptr_t A = MappingFactory::createInstance(SpaceX, SpaceY, matrix, false);
		Eigen::VectorXd vector(10);
		vector << 0.223629,-0.802647,0.072957,0.424295,0.047306,-0.303829,-0.832450,0.722259,-0.066522,-0.552865;
		const SpaceElement_ptr_t x =
				ElementCreator::create(SpaceX, vector);
		const SpaceElement_ptr_t compare = (*A)(x);
		//std::cout << compare << std::endl;

		Eigen::VectorXd expected(2);
		expected << 2.45722, 0.201275;
		const SpaceElement_ptr_t expected_vector =
				ElementCreator::create(SpaceY, expected);

		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
}
