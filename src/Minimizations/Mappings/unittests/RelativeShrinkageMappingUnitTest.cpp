/*
 * RelativeShrinkageMappingUnitTest.cpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "RelativeShrinkageMappingUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/Specifics/RelativeShrinkageMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( RelativeShrinkageMappingUnitTest );

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
	const double power = 2.;
	Eigen::VectorXd xtemp(10);
	xtemp << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929;
	{
		const double lambda = .9;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createRegularizedL1Instance(
						xtemp.innerSize(), lambda, power);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
		// Note that we have the duality mapping (the soft thresholding)
		// only in the dual space, not in the original back as this
		// space is not smooth and hence its duality mapping not single-
		// valued.
		const Mapping_ptr_t S = SpaceX->getDualSpace()->getDualityMapping();
		Eigen::VectorXd expected(10);
		expected << 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.;
		const SpaceElement_ptr_t compare = (*S)(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
	{
		const double lambda = .4;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createRegularizedL1Instance(
						xtemp.innerSize(), lambda, power);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
		const Mapping_ptr_t S = SpaceX->getDualSpace()->getDualityMapping();
		Eigen::VectorXd expected(10);
		expected <<
0.00000,-0.99878, 0.00000, 0.00000, 0.00000, 0.00000,-1.12907, 0.80146,
0.00000,-0.32482;
		const SpaceElement_ptr_t compare = (*S)(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
	{
		const double lambda = .1;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createRegularizedL1Instance(
						xtemp.innerSize(), lambda, power);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
		const Mapping_ptr_t S = SpaceX->getDualSpace()->getDualityMapping();
		Eigen::VectorXd expected(10);
		expected <<
1.04691,-6.99513, 0.00000, 2.64664, 0.00000,-1.72607,-7.51628, 6.20586,
0.00000,-4.29929;
		const SpaceElement_ptr_t compare = (*S)(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
	{
		const double lambda = .001;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createRegularizedL1Instance(
						xtemp.innerSize(), lambda, power);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
		const Mapping_ptr_t S = SpaceX->getDualSpace()->getDualityMapping();
		Eigen::VectorXd expected(10);
		expected <<
203.691,-798.513,55.042, 363.664,38.179,-271.607,-850.628, 719.586,
-57.074,-528.929;
		const SpaceElement_ptr_t compare = (*S)(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
//	{
//		const double lambda = .0;
//		const NormedSpace_ptr_t SpaceX =
//				NormedSpaceFactory::createRegularizedL1Instance(
//						xtemp.innerSize(), lambda, power);
//		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
//		const Mapping_ptr_t S = SpaceX->getDualSpace()->getDualityMapping();
//		Eigen::VectorXd expected(xtemp);
//		const SpaceElement_ptr_t compare = (*S)(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare << ".\n";
//		SpaceElement_ptr_t expected_vector =
//				ElementCreator::create(compare->getSpace(), expected);
//		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
//	}
}
