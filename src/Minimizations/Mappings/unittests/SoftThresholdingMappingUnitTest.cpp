/*
 * SoftThresholdingMappingUnitTest.cpp
 *
 *  Created on: Oct 13, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "SoftThresholdingMappingUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/Specifics/SoftThresholdingMapping.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SoftThresholdingMappingUnitTest );

// static entities
double SoftThresholdingMappingUnitTest::tolerance = 1e-4;


void SoftThresholdingMappingUnitTest::setUp()
{
}


void SoftThresholdingMappingUnitTest::tearDown()
{
}

/** We generate test vectors as follows via octave:
 *
 * x=ones(2,10)-2.*rand(2,10)
 * gval=zeros(2,10)
 * for i=1:10
 * 	gval(:,i)=SoftThresholdingMapping(x, p, power, 1e-6)
 * endfor
 * gval
 *
 *
 */

void SoftThresholdingMappingUnitTest::oneNorm()
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
		expected << 0.,-0.399513,0.,0.,0.,0.,-0.451628,0.320586,0.,-0.129929;
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
		expected << 0.104691,-0.699513,0.,0.264664,0.,-0.172607,-0.751628,0.620586,0.,-0.429929;
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
		expected << 0.203691,-0.798513,0.055042,0.363664,0.038179,-0.271607,-0.850628,0.719586,-0.057074,-0.528929;
		const SpaceElement_ptr_t compare = (*S)(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
	{
		const double lambda = .0;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createRegularizedL1Instance(
						xtemp.innerSize(), lambda, power);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX,xtemp);
		const Mapping_ptr_t S = SpaceX->getDualSpace()->getDualityMapping();
		Eigen::VectorXd expected(xtemp);
		const SpaceElement_ptr_t compare = (*S)(x);
//			std::cout << "Expecting " << expected.transpose()
//					<< " and got " << compare << ".\n";
		SpaceElement_ptr_t expected_vector =
				ElementCreator::create(compare->getSpace(), expected);
		CPPUNIT_ASSERT( expected_vector->isApprox(compare, tolerance)  );
	}
}
