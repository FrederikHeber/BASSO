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
 * RelativeShrinkageMappingUnitTest.cpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "RelativeShrinkageMappingUnitTest.hpp"

#include <boost/assign.hpp>

#include <Eigen/Dense>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/Specifics/RelativeShrinkageMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( RelativeShrinkageMappingUnitTest );

using namespace boost::assign;

// static entities
double RelativeShrinkageMappingUnitTest::tolerance = 1e-4;


void RelativeShrinkageMappingUnitTest::setUp()
{
}


void RelativeShrinkageMappingUnitTest::tearDown()
{
}

/** We generate test vectors as follows via octave:
 *
 * x=ones(2,10)-2.*rand(2,10)
 * gval=zeros(2,10)
 * for i=1:10
 * 	gval(:,i)=RelativeShrinkageMapping(x, p, power, 1e-6)
 * endfor
 * gval
 *
 *
 */

void RelativeShrinkageMappingUnitTest::oneNorm()
{
//	const double power = 2.;
	Eigen::VectorXd xtemp(10);
	xtemp << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929;
	{
		const double lambda = 1e+3;
		NormedSpaceFactory::args_t args;
		args += boost::any(lambda), boost::any(lambda);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						xtemp.innerSize(), "regularized_l1", args);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
		// Note that we have the duality mapping (the relative shrinkage)
		// only in the dual space, not in the original back as this
		// space is not smooth and hence its duality mapping not single-
		// valued.
		const RelativeShrinkageMapping &S = static_cast<RelativeShrinkageMapping&>
			(*SpaceX->getDualSpace()->getDualityMapping());
		const double coefficient = 0.00385832970297;
		const double compare_coefficient = S.getRelativeShrinkage(x);
		std::cout << "Expecting shrinkage coefficient " << coefficient
				<< " and got " << compare_coefficient << ".\n";
		CPPUNIT_ASSERT( fabs(coefficient - compare_coefficient) < BASSOTOLERANCE  );
		Eigen::VectorXd expected(10);
		expected <<
0.000200833,-0.000795655,5.21837e-05,0.000360806,3.53207e-05,-0.000268749,-0.00084777,0.000716728,-5.42157e-05,-0.000526071;
		const SpaceElement_ptr_t compare = S(x);
//		std::cout << "Expecting " << expected.transpose()
//				<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
	{
		const double lambda = 1e+1;
		NormedSpaceFactory::args_t args;
		args += boost::any(lambda), boost::any(lambda);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						xtemp.innerSize(), "regularized_l1", args);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
		const RelativeShrinkageMapping &S = static_cast<RelativeShrinkageMapping&>
			(*SpaceX->getDualSpace()->getDualityMapping());
		const double coefficient = 0.2211829375;
		const double compare_coefficient = S.getRelativeShrinkage(x);
//		std::cout << "Expecting shrinkage coefficient " << coefficient
//				<< " and got " << compare_coefficient << ".\n";
		CPPUNIT_ASSERT( fabs(coefficient - compare_coefficient) < BASSOTOLERANCE  );
		Eigen::VectorXd expected(10);
		expected <<
0,-0.057833,0,0.0143481,0,-0.00514241,-0.0630445,0.0499403,0,-0.0308746;
		const SpaceElement_ptr_t compare = S(x);
//		std::cout << "Expecting " << expected.transpose()
//				<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
	{
		const double lambda = .9;
		NormedSpaceFactory::args_t args;
		args += boost::any(lambda), boost::any(lambda);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						xtemp.innerSize(), "regularized_l1", args);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
		const RelativeShrinkageMapping &S = static_cast<RelativeShrinkageMapping&>
			(*SpaceX->getDualSpace()->getDualityMapping());
		const double coefficient = 0.6081351282051;
		const double compare_coefficient = S.getRelativeShrinkage(x);
//		std::cout << "Expecting shrinkage coefficient " << coefficient
//				<< " and got " << compare_coefficient << ".\n";
		CPPUNIT_ASSERT( fabs(coefficient - compare_coefficient) < BASSOTOLERANCE  );
		Eigen::VectorXd expected(10);
		expected <<
0,-0.2126420797721,0,0,0,0,-0.2705476353276,0.1249454131054,0,0.;
		const SpaceElement_ptr_t compare = S(x);
//		std::cout << "Expecting " << expected.transpose()
//				<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
	{
		const double lambda = .4;
		NormedSpaceFactory::args_t args;
		args += boost::any(lambda), boost::any(lambda);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						xtemp.innerSize(), "regularized_l1", args);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
		const RelativeShrinkageMapping &S = static_cast<RelativeShrinkageMapping&>
			(*SpaceX->getDualSpace()->getDualityMapping());
		const double coefficient = 0.6975667647059;
		const double compare_coefficient = S.getRelativeShrinkage(x);
//		std::cout << "Expecting shrinkage coefficient " << coefficient
//				<< " and got " << compare_coefficient << ".\n";
		CPPUNIT_ASSERT( fabs(coefficient - compare_coefficient) < BASSOTOLERANCE  );
		Eigen::VectorXd expected(10);
		expected <<
0,-0.254866,0,0,0,0,-0.385153,0.0575481,0,0.;
		const SpaceElement_ptr_t compare = S(x);
//		std::cout << "Expecting " << expected.transpose()
//				<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
	{
		const double lambda = .1;
		NormedSpaceFactory::args_t args;
		args += boost::any(lambda), boost::any(lambda);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						xtemp.innerSize(), "regularized_l1", args);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
		const RelativeShrinkageMapping &S = static_cast<RelativeShrinkageMapping&>
			(*SpaceX->getDualSpace()->getDualityMapping());
		const double coefficient = 0.7862576190476;
		const double compare_coefficient = S.getRelativeShrinkage(x);
//		std::cout << "Expecting shrinkage coefficient " << coefficient
//				<< " and got " << compare_coefficient << ".\n";
		CPPUNIT_ASSERT( fabs(coefficient - compare_coefficient) < BASSOTOLERANCE  );
		Eigen::VectorXd expected(10);
		expected <<
0,-0.1325538095238,0,0,0,0,-0.6537038095238,0,0,0.;
		const SpaceElement_ptr_t compare = S(x);
//		std::cout << "Expecting " << expected.transpose()
//				<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
	{
		const double lambda = .001;
		NormedSpaceFactory::args_t args;
		args += boost::any(lambda), boost::any(lambda);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						xtemp.innerSize(), "regularized_l1", args);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
		const RelativeShrinkageMapping &S = static_cast<RelativeShrinkageMapping&>
			(*SpaceX->getDualSpace()->getDualityMapping());
		const double coefficient = 0.8507772227772;
		const double compare_coefficient = S.getRelativeShrinkage(x);
//		std::cout << "Expecting shrinkage coefficient " << coefficient
//				<< " and got " << compare_coefficient << ".\n";
		CPPUNIT_ASSERT( fabs(coefficient - compare_coefficient) < BASSOTOLERANCE  );
		Eigen::VectorXd expected(10);
		expected <<
0,0,0,0,0,0,-0.8507772227772,0,0,0.;
		const SpaceElement_ptr_t compare = S(x);
//		std::cout << "Expecting " << expected.transpose()
//				<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
//	{
//		const double lambda = .0;
//		NormedSpaceFactory::args_t args;
//		args += boost::any(lambda), boost::any(lambda);
//		const NormedSpace_ptr_t SpaceX =
//				NormedSpaceFactory::create(
//						xtemp.innerSize(), "regularized_l1", args);
//		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
//		const RelativeShrinkageMapping &S = static_cast<RelativeShrinkageMapping&>
//			(*SpaceX->getDualSpace()->getDualityMapping());
//		const double coefficient = 0.6081351282051;
//		const double compare_coefficient = S.getRelativeShrinkage(x);
//		std::cout << std::setprecision(13) << "Expecting shrinkage coefficient " << coefficient
//				<< " and got " << compare_coefficient << ".\n";
//		CPPUNIT_ASSERT( fabs(coefficient - compare_coefficient) < BASSOTOLERANCE  );
//		Eigen::VectorXd expected(xtemp);
//		const SpaceElement_ptr_t compare = S(x);
//		std::cout << "Expecting " << expected.transpose()
//				<< " and got " << compare << ".\n";
//		SpaceElement_ptr_t expected_vector =
//				ElementCreator::create(compare->getSpace(), expected);
//		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
//	}
}
