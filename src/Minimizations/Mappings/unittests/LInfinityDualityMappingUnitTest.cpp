/*
 * LInfinityDualityMappingUnitTest.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LInfinityDualityMappingUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/Mappings/LInfinityDualityMapping.hpp"
#include "Minimizations/Norms/NormExceptions.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( LInfinityDualityMappingUnitTest );


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
	Eigen::MatrixXd X(2,10);
	X << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929,
			0.608996,0.567817,0.261505,-0.294843,-0.387682,-0.513624,-0.728372,-0.676635,-0.503763,-0.611381;
	{
		const double power = 1.;
		const LInfinityDualityMapping J_infty(power);
		Eigen::MatrixXd expected(2,10);
		expected << 0,-1,0,1,-0,-0,-1,1,-0,-0,
				1,0,1,0,-1,-1,0,0,-1,-1;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_infty(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 1.1;
		const LInfinityDualityMapping J_infty(power);
		Eigen::MatrixXd expected(2,10);
		expected << 0.00000,-0.97787,0.00000,0.90404,-0.00000,-0.00000,-0.98407,0.96776,-0.00000,-0.00000,
				0.95162,0.00000,0.87448,0.00000,-0.90959,-0.93554,0.00000,0.00000,-0.93373,-0.95199;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_infty(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 2.;
		const LInfinityDualityMapping J_infty(power);
		Eigen::MatrixXd expected(2,10);
		expected << 0.00000,-0.79951,0.00000,0.36466,-0.00000,-0.00000,-0.85163,0.72059,-0.00000,-0.00000,
				0.60900,0.00000,0.26151,0.00000,-0.38768,-0.51362,0.00000,0.00000,-0.50376,-0.61138;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_infty(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 10.;
		const LInfinityDualityMapping J_infty(power);
		Eigen::MatrixXd expected(2,10);
		expected << 0.00000,-0.13348,0.00000,0.00011,-0.00000,-0.00000,-0.23564,0.05238,-0.00000,-0.00000,
				0.01152,0.00000,0.00001,0.00000,-0.00020,-0.00249,0.00000,0.00000,-0.00209,-0.01193;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_infty(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( (expected.col(i) - compare).lpNorm<2>() < 1e-4  );
		}
	}
}

void LInfinityDualityMappingUnitTest::setTolerance()
{
	LInfinityDualityMapping J_infty(2.);
	CPPUNIT_ASSERT_EQUAL(BASSOTOLERANCE, J_infty.tolerance);
	const double value = 1e-1;
	J_infty.setTolerance(value);
	CPPUNIT_ASSERT_EQUAL(value, J_infty.tolerance);
}
