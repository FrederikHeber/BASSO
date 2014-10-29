/*
 * L1DualityMappingUnitTest.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "L1DualityMappingUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/Mappings/L1DualityMapping.hpp"
#include "Minimizations/Norms/NormExceptions.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( L1DualityMappingUnitTest );


void L1DualityMappingUnitTest::setUp()
{
}


void L1DualityMappingUnitTest::tearDown()
{
}

/** We generate test vectors as follows via octave:
 *
 * x=ones(2,10)-2.*rand(2,10)
 * gval=zeros(2,10)
 * for i=1:10
 * 	gval(:,i)=L1DualityMapping(x, p, power, 1e-6)
 * endfor
 * gval
 *
 *
 */

void L1DualityMappingUnitTest::oneNorm()
{
	const double p = 1.;
	Eigen::MatrixXd X(2,10);
	X << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929,
			0.608996,0.567817,0.261505,-0.294843,-0.387682,-0.513624,-0.728372,-0.676635,-0.503763,-0.611381;
	{
		const double power = .9;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						X.innerSize(), p, power);
		const Mapping &J_1 = *SpaceX->getSpace()->getDualityMapping();
		Eigen::MatrixXd expected(2,10);
		expected << 1.02083,-0.96920,1.12155,1.04250,1.08886,-1.02434,-0.95529,0.96710,-1.05935,-0.98687,
				1.02083,0.96920,1.12155,-1.04250,-1.08886,-1.02434,-0.95529,-0.96710,-1.05935,-0.98687;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_1(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 1.;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						X.innerSize(), p, power);
		const Mapping &J_1 = *SpaceX->getSpace()->getDualityMapping();
		Eigen::MatrixXd expected(2,10);
		expected << 1,-1,1,1,1,-1,-1,1,-1,-1,
				1,1,1,-1,-1,-1,-1,-1,-1,-1;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_1(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 1.1;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						X.innerSize(), p, power);
		const Mapping &J_1 = *SpaceX->getSpace()->getDualityMapping();
		Eigen::MatrixXd expected(2,10);
		expected << 0.97959,-1.03178,0.89162,0.95923,0.91839,-0.97624,-1.04680,1.03401,-0.94398,-1.01331,
				0.97959,1.03178,0.89162,-0.95923,-0.91839,-0.97624,-1.04680,-1.03401,-0.94398,-1.01331;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_1(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 2.;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						X.innerSize(), p, power);
		const Mapping &J_1 = *SpaceX->getSpace()->getDualityMapping();
		Eigen::MatrixXd expected(2,10);
		expected << 0.81369,-1.36733,0.31755,0.65951,0.42686,-0.78623,-1.58000,1.39722,-0.56184,-1.14131,
				0.81369,1.36733,0.31755,-0.65951,-0.42686,-0.78623,-1.58000,-1.39722,-0.56184,-1.14131;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_1(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 10.;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						X.innerSize(), p, power);
		const Mapping &J_1 = *SpaceX->getSpace()->getDualityMapping();
		Eigen::MatrixXd expected(2,10);
		expected << 1.5636e-01,-1.6706e+01,3.2830e-05,2.3603e-02,4.7052e-04,-1.1480e-01,-6.1364e+01,2.0295e+01,-5.5782e-03,-3.2857e+00,
				1.5636e-01,1.6706e+01,3.2830e-05,-2.3603e-02,-4.7052e-04,-1.1480e-01,-6.1364e+01,-2.0295e+01,-5.5782e-03,-3.2857e+00;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_1(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
}

void L1DualityMappingUnitTest::setTolerance()
{
	const double p = 1.;
	const double power = 2.;
	const NormedSpace_ptr_t SpaceX =
			NormedSpaceFactory::createLpInstance(
					10, p, power);
	const PowerTypeDualityMapping &J_1 =
			dynamic_cast<const PowerTypeDualityMapping &>(
					*SpaceX->getSpace()->getDualityMapping());
	CPPUNIT_ASSERT_EQUAL(BASSOTOLERANCE, J_1.getTolerance());
	const double value = 1e-1;
	J_1.setTolerance(value);
	CPPUNIT_ASSERT_EQUAL(value, J_1.getTolerance());
}
