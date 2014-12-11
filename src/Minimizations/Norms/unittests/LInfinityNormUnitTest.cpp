/*
 * LInfinityNormUnitTest.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#include "LInfinityNormUnitTest.hpp"

#include <Eigen/Dense>
#include <limits>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Norms/LInfinityNorm.hpp"
#include "Minimizations/Norms/NormExceptions.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( LInfinityNormUnitTest );


void LInfinityNormUnitTest::setUp()
{
}


void LInfinityNormUnitTest::tearDown()
{
}

void LInfinityNormUnitTest::inftyNorm()
{
	const double p = std::numeric_limits<double>::infinity();
	{
		Eigen::MatrixXd X(2,10);
		X << 0.038619,-0.888852,0.206035,-0.192165,0.157826,-0.129335,-0.057939,-0.735626,-0.205108,-0.562922,
				0.520878,0.751919,-0.780070,-0.446509,0.096160,-0.903741,0.880418,-0.997373,0.345203,0.204716;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						X.innerSize(), p, 2.);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 0.52088,0.88885,0.78007,0.44651,0.15783,0.90374,0.88042,0.99737,0.34520,0.56292;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < 1e-4  );
		}
	}
	{
		Eigen::MatrixXd X(5,10);
		X << 0.2465831,-0.8789274,0.8060486,-0.2690611,0.4783144,-0.0882388,0.5841966,-0.9873399,0.1025274,0.1147468,
			0.7420882,-0.3929283,-0.8972077,-0.0136248,0.9230047,0.7008641,0.3326720,-0.6046968,-0.0304502,0.5890823,
			0.5836364,0.7425789,0.1059589,0.3170492,0.6488256,-0.6131489,0.3917273,-0.4766553,0.4321104,-0.2287513,
			0.0010107,0.8335353,0.7493321,0.1686532,-0.9022343,-0.1451915,-0.2955828,-0.8714578,-0.6201104,-0.4218609,
			0.9747178,-0.9666331,0.2735608,0.1137258,0.1620735,-0.4401804,-0.0051763,-0.6895868,-0.2982194,0.9003390;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						X.innerSize(), p, 2.);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 0.97472,0.96663,0.89721,0.31705,0.92300,0.70086,0.58420,0.98734,0.62011,0.90034;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < 1e-4  );
		}
	}
	{
		Eigen::MatrixXd X(10,10);
		X << 4.1486e-01,9.6560e-01,-9.2229e-01,-1.6217e-01,4.4662e-01,-4.7275e-01,7.2439e-01,2.9047e-01,5.6243e-01,-4.1897e-02,
			-7.1267e-01,-9.4491e-01,4.8262e-01,-5.7532e-01,5.8245e-02,-6.3297e-01,7.8772e-02,-3.0798e-02,-9.1078e-01,6.7228e-01,
			9.8346e-01,3.6610e-01,-6.2597e-02,-8.4511e-04,-8.6594e-01,-2.7087e-01,3.4065e-01,4.2926e-02,8.3941e-01,-7.1284e-01,
			2.1323e-01,-8.8896e-01,6.9595e-01,-6.6868e-02,-2.2635e-01,-1.1158e-01,4.8951e-01,7.4771e-01,6.0980e-01,5.6190e-01,
			-7.5020e-01,4.6583e-01,3.1282e-01,-1.6002e-02,-1.6480e-01,7.0506e-02,2.4843e-01,7.2700e-01,9.8683e-01,5.9212e-02,
			8.8330e-01,8.4138e-01,-3.1399e-01,-8.5381e-02,3.4769e-01,-1.4253e-01,1.9348e-01,3.6669e-01,-5.2582e-01,-5.1170e-01,
			7.5041e-04,-3.9968e-01,3.8499e-01,5.9907e-01,-9.5389e-01,-8.3976e-01,6.2231e-01,-1.3057e-01,-1.8735e-01,9.6702e-01,
			-5.3199e-01,-1.8225e-01,-3.3247e-01,8.7822e-01,-3.9038e-01,7.2615e-01,-8.6705e-01,2.9549e-02,1.8063e-02,-4.1061e-01,
			-5.1152e-01,-8.3150e-01,2.7693e-01,-9.7872e-01,9.7519e-01,-1.7073e-01,7.9553e-01,-4.0667e-01,-9.7718e-01,8.2855e-01,
			1.4308e-01,8.0489e-02,-5.3567e-01,-3.6793e-02,7.2730e-01,1.9200e-01,-4.9364e-01,6.0502e-01,2.9298e-01,-5.2860e-01;
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						X.innerSize(), p, 2.);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 0.98346,0.96560,0.92229,0.97872,0.97519,0.83976,0.86705,0.74771,0.98683,0.96702;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < 1e-4  );
		}
	}
}
