/*
 * LpNormUnitTest.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "LpNormUnitTest.hpp"

#include <boost/assign.hpp>

#include <Eigen/Dense>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/NormExceptions.hpp"
#include "Minimizations/Norms/NormFactory.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( LpNormUnitTest );

using namespace boost::assign;

void LpNormUnitTest::setUp()
{
}


void LpNormUnitTest::tearDown()
{
}

/** We generate test vectors as follows via octave:
 *
 * x=ones(2,10)-2.*rand(2,10)
 * norm(x, p, "cols")
 *
 */



void LpNormUnitTest::otherNorm()
{
	Eigen::MatrixXd X(2,10);
	X << 0.038619,-0.888852,0.206035,-0.192165,0.157826,-0.129335,-0.057939,-0.735626,-0.205108,-0.562922,
			0.520878,0.751919,-0.780070,-0.446509,0.096160,-0.903741,0.880418,-0.997373,0.345203,0.204716;
	{
		const double p = 1.1;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 0.54787,1.54111,0.94243,0.60454,0.23918,1.00005,0.92045,1.62903,0.51839,0.72887;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < 1e-4  );
		}
	}
	{
		const double p = 1.5;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 0.52787,1.30455,0.84915,0.52703,0.20456,0.93607,0.89030,1.38333,0.44386,0.64248;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < 1e-4  );
		}
	}
	{
		const double p = 5.64763763;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 0.52088,0.94207,0.78015,0.44718,0.15949,0.90374,0.88042,1.02691,0.34837,0.56325;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < 1e-4  );
		}
	}
}

void LpNormUnitTest::twoNorm()
{
	const double p = 2.;
	{
		Eigen::MatrixXd X(2,10);
		X << 0.038619,-0.888852,0.206035,-0.192165,0.157826,-0.129335,-0.057939,-0.735626,-0.205108,-0.562922,
				0.520878,0.751919,-0.780070,-0.446509,0.096160,-0.903741,0.880418,-0.997373,0.345203,0.204716;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 0.52231,1.16423,0.80682,0.48610,0.18481,0.91295,0.88232,1.23931,0.40154,0.59899;
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
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 1.37920,1.76281,1.44992,0.46312,1.53036,1.04393,0.83235,1.67412,0.81953,1.18368;
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
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 1.8976,2.1367,1.5462,1.5680,1.9220,1.4261,1.7324,1.3633,2.1328,1.9007;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < 1e-4  );
		}
	}
}

void LpNormUnitTest::threeNorm()
{
	const double p = 3.;
	{
		Eigen::MatrixXd X(2,10);
		X << 0.038619,-0.888852,0.206035,-0.192165,0.157826,-0.129335,-0.057939,-0.735626,-0.205108,-0.562922,
				0.520878,0.751919,-0.780070,-0.446509,0.096160,-0.903741,0.880418,-0.997373,0.345203,0.204716;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 0.52095,1.04077,0.78483,0.45807,0.16893,0.90462,0.88050,1.11608,0.36782,0.57181;
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
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 1.15693,1.38059,1.19075,0.38624,1.24021,0.87234,0.68551,1.31649,0.70255,1.00754;
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
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 1.4083,1.5656,1.1609,1.2653,1.4480,1.1109,1.2758,1.0529,1.5636,1.3869;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < 1e-4  );
		}
	}
}

void LpNormUnitTest::fourNorm()
{
	const double p = 4.;
	{
		Eigen::MatrixXd X(2,10);
		X << 0.038619,-0.888852,0.206035,-0.192165,0.157826,-0.129335,-0.057939,-0.735626,-0.205108,-0.562922,
				0.520878,0.751919,-0.780070,-0.446509,0.096160,-0.903741,0.880418,-0.997373,0.345203,0.204716;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 0.52088,0.98566,0.78102,0.45029,0.16300,0.90384,0.88042,1.06415,0.35549,0.56537;
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
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 1.07302,1.22887,1.08603,0.35743,1.12795,0.80535,0.63236,1.17996,0.66087,0.94930;
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
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 1.23144,1.35409,1.03648,1.15041,1.27653,0.99849,1.11154,0.93939,1.35471,1.20136;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < 1e-4  );
		}
	}
}

void LpNormUnitTest::elevenNorm()
{
	const double p = 11.;
	{
		Eigen::MatrixXd X(2,10);
		X << 0.038619,-0.888852,0.206035,-0.192165,0.157826,-0.129335,-0.057939,-0.735626,-0.205108,-0.562922,
				0.520878,0.751919,-0.780070,-0.446509,0.096160,-0.903741,0.880418,-0.997373,0.345203,0.204716;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 0.52088,0.90084,0.78007,0.44651,0.15789,0.90374,0.88042,1.00051,0.34531,0.56292;
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
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 0.97934,1.00896,0.92778,0.32149,0.97366,0.71448,0.58498,1.00956,0.62118,0.90113;
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
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(2.);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(
						X.innerSize(), "lp", args);
		const Norm &norm = *SpaceX->getNorm();
		Eigen::VectorXd expected(10);
		expected << 1.01326,1.05338,0.92628,1.00312,1.04300,0.85687,0.90262,0.79006,1.07209,0.98618;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < 1e-4  );
		}
	}
}
