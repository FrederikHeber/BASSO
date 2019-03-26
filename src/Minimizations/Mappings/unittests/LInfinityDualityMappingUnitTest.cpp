/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
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
 * LInfinityDualityMappingUnitTest.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LInfinityDualityMappingUnitTest.hpp"

#include <boost/assign.hpp>

#include <Eigen/Dense>
#include <limits>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/Specifics/LInfinityDualityMapping.hpp"
#include "Minimizations/Norms/NormExceptions.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( LInfinityDualityMappingUnitTest );

using namespace boost::assign;

// static entities
double LInfinityDualityMappingUnitTest::tolerance = 1e-4;


void LInfinityDualityMappingUnitTest::setUp()
{
}


void LInfinityDualityMappingUnitTest::tearDown()
{
}

/** We generate test vectors as follows via octave:
 *
 * x=ones(2,10)-2.*rand(2,10)
 * gval=zeros(2,10)
 * for i=1:10
 * 	gval(:,i)=LInfinityDualityMapping(x, p, power, 1e-6)
 * endfor
 * gval
 *
 *
 */

void LInfinityDualityMappingUnitTest::inftyNorm()
{
	const double p = std::numeric_limits<double>::infinity();
	Eigen::MatrixXd X(2,10);
	X << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929,
			0.608996,0.567817,0.261505,-0.294843,-0.387682,-0.513624,-0.728372,-0.676635,-0.503763,-0.611381;
	{
		const double power = 1.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Mapping &J_infty = *SpaceX->getSpace()->getDualityMapping();
		Eigen::MatrixXd expected(2,10);
		expected << 0,-1,0,1,-0,-0,-1,1,-0,-0,
				1,0,1,0,-1,-1,0,0,-1,-1;
		for (size_t i=0; i<10; ++i) {
			const SpaceElement_ptr_t x =
					ElementCreator::create(SpaceX, X.col(i));
			const SpaceElement_ptr_t compare = J_infty(x);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare << ".\n";
			SpaceElement_ptr_t expected_vector =
					ElementCreator::create(compare->getSpace(), expected.col(i));
			CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
		}
	}
	{
		const double power = 1.1;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Mapping &J_infty = *SpaceX->getSpace()->getDualityMapping();
		Eigen::MatrixXd expected(2,10);
		expected << 0.00000,-0.97787,0.00000,0.90404,-0.00000,-0.00000,-0.98407,0.96776,-0.00000,-0.00000,
				0.95162,0.00000,0.87448,0.00000,-0.90959,-0.93554,0.00000,0.00000,-0.93373,-0.95199;
		for (size_t i=0; i<10; ++i) {
			const SpaceElement_ptr_t x =
					ElementCreator::create(SpaceX, X.col(i));
			const SpaceElement_ptr_t compare = J_infty(x);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare << ".\n";
			SpaceElement_ptr_t expected_vector =
					ElementCreator::create(compare->getSpace(), expected.col(i));
			CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
		}
	}
	{
		const double power = 2.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Mapping &J_infty = *SpaceX->getSpace()->getDualityMapping();
		Eigen::MatrixXd expected(2,10);
		expected << 0.00000,-0.79951,0.00000,0.36466,-0.00000,-0.00000,-0.85163,0.72059,-0.00000,-0.00000,
				0.60900,0.00000,0.26151,0.00000,-0.38768,-0.51362,0.00000,0.00000,-0.50376,-0.61138;
		SpaceElement_ptr_t x = SpaceX->createElement();
		for (size_t i=0; i<10; ++i) {
			const SpaceElement_ptr_t x =
					ElementCreator::create(SpaceX, X.col(i));
			const SpaceElement_ptr_t compare = J_infty(x);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare << ".\n";
			SpaceElement_ptr_t expected_vector =
					ElementCreator::create(compare->getSpace(), expected.col(i));
			CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
		}
	}
	{
		const double power = 10.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Mapping &J_infty = *SpaceX->getSpace()->getDualityMapping();
		Eigen::MatrixXd expected(2,10);
		expected << 0.00000,-0.13348,0.00000,0.00011,-0.00000,-0.00000,-0.23564,0.05238,-0.00000,-0.00000,
				0.01152,0.00000,0.00001,0.00000,-0.00020,-0.00249,0.00000,0.00000,-0.00209,-0.01193;
		SpaceElement_ptr_t x = SpaceX->createElement();
		for (size_t i=0; i<10; ++i) {
			const SpaceElement_ptr_t x =
					ElementCreator::create(SpaceX, X.col(i));
			const SpaceElement_ptr_t compare = J_infty(x);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare << ".\n";
			SpaceElement_ptr_t expected_vector =
					ElementCreator::create(compare->getSpace(), expected.col(i));
			// compare with absolute precision due to components less than 1
			CPPUNIT_ASSERT( (expected_vector-compare)->isZero(tolerance)  );
		}
	}
}

void LInfinityDualityMappingUnitTest::setTolerance()
{
	const double p = std::numeric_limits<double>::infinity();
	const double power = 2.;
	NormedSpaceFactory::args_t args;
	args += boost::any(p), boost::any(power);
	const NormedSpace_ptr_t SpaceX =
			NormedSpaceFactory::create(
					10, "lp", args);
	const PowerTypeDualityMapping &J_infty =
			dynamic_cast<const PowerTypeDualityMapping &>(
					*SpaceX->getSpace()->getDualityMapping());
	Eigen::MatrixXd expected(2,10);
	CPPUNIT_ASSERT_EQUAL(BASSOTOLERANCE, J_infty.getTolerance());
	const double value = 1e-1;
	J_infty.setTolerance(value);
	CPPUNIT_ASSERT_EQUAL(value, J_infty.getTolerance());
}
