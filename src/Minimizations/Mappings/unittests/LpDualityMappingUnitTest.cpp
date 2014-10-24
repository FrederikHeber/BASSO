/*
 * LpDualityMappingUnitTest.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "LpDualityMappingUnitTest.hpp"

#include <Eigen/Dense>

#include "Minimizations/Mappings/LpDualityMapping.hpp"
#include "Minimizations/Norms/NormExceptions.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( LpDualityMappingUnitTest );


void LpDualityMappingUnitTest::setUp()
{
}


void LpDualityMappingUnitTest::tearDown()
{
}

void LpDualityMappingUnitTest::throwTest()
{
	// we check that assertion is thrown for invalid p value
//	std::cout << "The following assertion is intended and does not indicate a failure of the test." << std::endl;
	CPPUNIT_ASSERT_THROW(
			LpDualityMapping J_illegal(-0.5),
			NormIllegalValue_exception );
}

/** We generate test vectors as follows via octave:
 *
 * x=ones(2,10)-2.*rand(2,10)
 * gval=zeros(2,10)
 * for i=1:10
 * 	gval(:,i)=LpDualityMapping(x, p, power, 1e-6)
 * endfor
 * gval
 *
 *
 */

void LpDualityMappingUnitTest::oneNorm()
{
	const double p = 1.;
	LpDualityMapping J_p(p);
	Eigen::MatrixXd X(2,10);
	X << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929,
			0.608996,0.567817,0.261505,-0.294843,-0.387682,-0.513624,-0.728372,-0.676635,-0.503763,-0.611381;
	{
		const double power = .9;
		Eigen::MatrixXd expected(2,10);
		expected << 1.02083,-0.96920,1.12155,1.04250,1.08886,-1.02434,-0.95529,0.96710,-1.05935,-0.98687,
				1.02083,0.96920,1.12155,-1.04250,-1.08886,-1.02434,-0.95529,-0.96710,-1.05935,-0.98687;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 1.;
		Eigen::MatrixXd expected(2,10);
		expected << 1,-1,1,1,1,-1,-1,1,-1,-1,
				1,1,1,-1,-1,-1,-1,-1,-1,-1;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 1.1;
		Eigen::MatrixXd expected(2,10);
		expected << 0.97959,-1.03178,0.89162,0.95923,0.91839,-0.97624,-1.04680,1.03401,-0.94398,-1.01331,
				0.97959,1.03178,0.89162,-0.95923,-0.91839,-0.97624,-1.04680,-1.03401,-0.94398,-1.01331;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 2.;
		Eigen::MatrixXd expected(2,10);
		expected << 0.81369,-1.36733,0.31755,0.65951,0.42686,-0.78623,-1.58000,1.39722,-0.56184,-1.14131,
				0.81369,1.36733,0.31755,-0.65951,-0.42686,-0.78623,-1.58000,-1.39722,-0.56184,-1.14131;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 10.;
		Eigen::MatrixXd expected(2,10);
		expected << 1.5636e-01,-1.6706e+01,3.2830e-05,2.3603e-02,4.7052e-04,-1.1480e-01,-6.1364e+01,2.0295e+01,-5.5782e-03,-3.2857e+00,
				1.5636e-01,1.6706e+01,3.2830e-05,-2.3603e-02,-4.7052e-04,-1.1480e-01,-6.1364e+01,-2.0295e+01,-5.5782e-03,-3.2857e+00;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
}

void LpDualityMappingUnitTest::twoNorm()
{
	const double p = 2.;
	LpDualityMapping J_p(p);
	Eigen::MatrixXd X(2,10);
	X << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929,
			0.608996,0.567817,0.261505,-0.294843,-0.387682,-0.513624,-0.728372,-0.676635,-0.503763,-0.611381;
	{
		const double power = 1.;
		Eigen::MatrixXd expected(2,10);
		expected << 0.31860,-0.81530,0.20955,0.77762,0.10055,-0.46881,-0.75996,0.72899,-0.11452,-0.65498,
				0.94789,0.57903,0.97780,-0.62873,-0.99493,-0.88330,-0.64997,-0.68453,-0.99342,-0.75565;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 1.8;
		Eigen::MatrixXd expected(2,10);
		expected << 0.223629,-0.802647,0.072957,0.424295,0.047306,-0.303829,-0.832450,0.722259,-0.066522,-0.552865,
				0.665339,0.570043,0.340435,-0.343057,-0.468101,-0.572450,-0.711969,-0.678206,-0.577042,-0.637843;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 2.;
		Eigen::MatrixXd expected(2,10);
		expected << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929,
				0.608996,0.567817,0.261505,-0.294843,-0.387682,-0.513624,-0.728372,-0.676635,-0.503763,-0.611381;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 2.1;
		Eigen::MatrixXd expected(2,10);
		expected << 0.195832,-0.797951,0.049118,0.338069,0.035655,-0.258221,-0.861382,0.719751,-0.054261,-0.518820,
				0.582640,0.566707,0.229194,-0.273340,-0.352812,-0.486518,-0.736714,-0.675851,-0.470691,-0.598565;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 10.;
		Eigen::MatrixXd expected(2,10);
		expected << 5.9422e-03,-6.8371e-01,1.4668e-06,8.5288e-04,2.0821e-05,-3.5632e-03,-2.1180e+00,6.5676e-01,-2.5394e-04,-9.7309e-02,
				1.7679e-02,4.8557e-01,6.8442e-06,-6.8959e-04,-2.0603e-04,-6.7135e-03,-1.8115e+00,-6.1670e-01,-2.2028e-03,-1.1227e-01;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
}

void LpDualityMappingUnitTest::fourNorm()
{
	const double p = 4.;
	LpDualityMapping J_p(p);
	Eigen::MatrixXd X(2,10);
	X << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929,
			0.608996,0.567817,0.261505,-0.294843,-0.387682,-0.513624,-0.728372,-0.676635,-0.503763,-0.611381;
	{
		const double power = 1.;
		Eigen::MatrixXd expected(2,10);
		expected << 0.0376116,-0.8436665,0.0098269,0.7657735,0.0010320,-0.1411895,-0.7251088,0.6496079,-0.0015318,-0.4655289,
				0.9905337,0.3022173,0.9984210,-0.4047568,-0.9999218,-0.9443374,-0.4536418,-0.5378451,-0.9998676,-0.7148738;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 2.;
		Eigen::MatrixXd expected(2,10);
		expected << 2.2978e-02,-7.1385e-01,2.5711e-03,3.0523e-01,4.0012e-04,-7.3916e-02,-6.8736e-01,5.4049e-01,-7.7171e-04,-3.1831e-01,
				6.0515e-01,2.5571e-01,2.6123e-01,-1.6133e-01,-3.8766e-01,-4.9438e-01,-4.3003e-01,-4.4750e-01,-5.0372e-01,-4.8880e-01;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 3.5;
		Eigen::MatrixXd expected(2,10);
		expected << 1.0972e-02,-5.5560e-01,3.4410e-04,7.6810e-02,9.6587e-05,-2.7999e-02,-6.3439e-01,4.1019e-01,-2.7594e-04,-1.7997e-01,
				2.8897e-01,1.9902e-01,3.4961e-02,-4.0598e-02,-9.3580e-02,-1.8727e-01,-3.9689e-01,-3.3962e-01,-1.8012e-01,-2.7637e-01;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 4.;
		Eigen::MatrixXd expected(2,10);
		expected << 8.5762e-03,-5.1107e-01,1.7601e-04,4.8493e-02,6.0140e-05,-2.0259e-02,-6.1766e-01,3.7416e-01,-1.9586e-04,-1.4882e-01,
				2.2586e-01,1.8307e-01,1.7883e-02,-2.5631e-02,-5.8268e-02,-1.3550e-01,-3.8642e-01,-3.0979e-01,-1.2784e-01,-2.2853e-01;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 4.2;
		Eigen::MatrixXd expected(2,10);
		expected << 7.7713e-03,-4.9427e-01,1.3461e-04,4.0344e-02,4.9757e-05,-1.7799e-02,-6.1109e-01,3.6065e-01,-1.7076e-04,-1.3792e-01,
				2.0466e-01,1.7706e-01,1.3677e-02,-2.1324e-02,-4.8209e-02,-1.1905e-01,-3.8231e-01,-2.9860e-01,-1.1146e-01,-2.1180e-01;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 20.;
		Eigen::MatrixXd expected(2,10);
		expected << 3.2297e-06,-3.5273e-02,8.4895e-14,1.9683e-08,1.5666e-11,-6.4505e-07,-2.6258e-01,1.9734e-02,-3.3719e-09,-3.3969e-04,
				8.5056e-05,1.2635e-02,8.6254e-12,-1.0404e-08,-1.5178e-08,-4.3144e-06,-1.6428e-01,-1.6339e-02,-2.2009e-06,-5.2164e-04;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
}

void LpDualityMappingUnitTest::elevenNorm()
{
	const double p = 11.;
	LpDualityMapping J_p(p);
	Eigen::MatrixXd X(2,10);
	X << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929,
			0.608996,0.567817,0.261505,-0.294843,-0.387682,-0.513624,-0.728372,-0.676635,-0.503763,-0.611381;
	{
		const double power = 1.;
		Eigen::MatrixXd expected(2,10);
		expected << 1.8401e-05,-9.7938e-01,2.0433e-07,9.1964e-01,1.1112e-10,-1.7723e-03,-8.6089e-01,6.9151e-01,-4.1453e-10,-2.0166e-01,
				9.9999e-01,3.1973e-02,1.0000e+00,-1.0980e-01,-1.0000e+00,-9.9914e-01,-1.8029e-01,-3.6854e-01,-1.0000e+00,-8.4249e-01;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 2.;
		Eigen::MatrixXd expected(2,10);
		expected << 1.1206e-05,-7.8466e-01,5.3434e-08,3.3818e-01,4.3078e-11,-9.1039e-04,-7.4422e-01,5.1702e-01,-2.0882e-10,-1.2542e-01,
				6.0899e-01,2.5616e-02,2.6150e-01,-4.0376e-02,-3.8768e-01,-5.1323e-01,-1.5586e-01,-2.7554e-01,-5.0376e-01,-5.2399e-01;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 10.9;
		Eigen::MatrixXd expected(2,10);
		expected << 1.3568e-07,-1.0911e-01,3.4945e-13,4.5960e-05,9.3688e-15,-2.4227e-06,-2.0362e-01,3.8859e-02,-4.6730e-13,-1.8315e-03,
				7.3737e-03,3.5621e-03,1.7102e-06,-5.4873e-06,-8.4316e-05,-1.3658e-03,-4.2644e-02,-2.0710e-02,-1.1273e-03,-7.6515e-03;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 11.;
		Eigen::MatrixXd expected(2,10);
		expected << 1.2912e-07,-1.0672e-01,3.0559e-13,4.1584e-05,8.5218e-15,-2.2666e-06,-2.0068e-01,3.7745e-02,-4.3633e-13,-1.7465e-03,
				7.0169e-03,3.4840e-03,1.4955e-06,-4.9649e-06,-7.6693e-05,-1.2778e-03,-4.2027e-02,-2.0116e-02,-1.0526e-03,-7.2966e-03;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 11.5;
		Eigen::MatrixXd expected(2,10);
		expected << 1.0076e-07,-9.5526e-02,1.5627e-13,2.5217e-05,5.3060e-15,-1.6245e-06,-1.8659e-01,3.2637e-02,-3.0969e-13,-1.3774e-03,
				5.4759e-03,3.1185e-03,7.6479e-07,-3.0107e-06,-4.7752e-05,-9.1579e-04,-3.9076e-02,-1.7394e-02,-7.4709e-04,-5.7544e-03;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
}

void LpDualityMappingUnitTest::inftyNorm()
{
	const double p = LpNorm::Infinity;
	LpDualityMapping J_p(p);
	Eigen::MatrixXd X(2,10);
	X << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929,
			0.608996,0.567817,0.261505,-0.294843,-0.387682,-0.513624,-0.728372,-0.676635,-0.503763,-0.611381;
	{
		const double power = 1.;
		Eigen::MatrixXd expected(2,10);
		expected << 0,-1,0,1,-0,-0,-1,1,-0,-0,
				1,0,1,0,-1,-1,0,0,-1,-1;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 1.1;
		Eigen::MatrixXd expected(2,10);
		expected << 0.00000,-0.97787,0.00000,0.90404,-0.00000,-0.00000,-0.98407,0.96776,-0.00000,-0.00000,
				0.95162,0.00000,0.87448,0.00000,-0.90959,-0.93554,0.00000,0.00000,-0.93373,-0.95199;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 2.;
		Eigen::MatrixXd expected(2,10);
		expected << 0.00000,-0.79951,0.00000,0.36466,-0.00000,-0.00000,-0.85163,0.72059,-0.00000,-0.00000,
				0.60900,0.00000,0.26151,0.00000,-0.38768,-0.51362,0.00000,0.00000,-0.50376,-0.61138;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double power = 10.;
		Eigen::MatrixXd expected(2,10);
		expected << 0.00000,-0.13348,0.00000,0.00011,-0.00000,-0.00000,-0.23564,0.05238,-0.00000,-0.00000,
				0.01152,0.00000,0.00001,0.00000,-0.00020,-0.00249,0.00000,0.00000,-0.00209,-0.01193;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( (expected.col(i) - compare).lpNorm<2>() < 1e-4  );
		}
	}
}

void LpDualityMappingUnitTest::otherNorm()
{
	Eigen::MatrixXd X(2,10);
	X << 0.204691,-0.799513,0.056042,0.364664,0.039179,-0.272607,-0.851628,0.720586,-0.058074,-0.529929,
			0.608996,0.567817,0.261505,-0.294843,-0.387682,-0.513624,-0.728372,-0.676635,-0.503763,-0.611381;
	{
		const double p = 1.1;
		const double power = p/(p-1.);
		LpDualityMapping J_p(p);
		Eigen::MatrixXd expected(2,10);
		expected << 6.7370e-02,-1.1768e+01,5.8495e-06,7.9058e-03,1.2227e-04,-4.5604e-02,-4.8991e+01,1.4229e+01,-1.8862e-03,-1.8659e+00,
				7.5131e-02,1.1372e+01,6.8236e-06,-7.7396e-03,-1.5377e-04,-4.8586e-02,-4.8231e+01,-1.4140e+01,-2.3411e-03,-1.8928e+00;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double p = 1.5;
		const double power = p/(p-1.);
		LpDualityMapping J_p(p);
		Eigen::MatrixXd expected(2,10);
		expected << 0.256915,-1.021804,0.034798,0.229659,0.049314,-0.266507,-1.298931,0.991715,-0.089537,-0.628823,
				0.443146,0.861111,0.075169,-0.206506,-0.155126,-0.365816,-1.201262,-0.960995,-0.263710,-0.675422;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double p = 5.;
		const double power = p/(p-1.);
		LpDualityMapping J_p(p);
		Eigen::MatrixXd expected(2,10);
		expected << 1.1238e-02,-8.3485e-01,1.5078e-03,6.2202e-01,8.2305e-05,-6.5131e-02,-7.2415e-01,6.1077e-01,-1.4879e-04,-3.7023e-01,
				8.8056e-01,2.1239e-01,7.1486e-01,-2.6582e-01,-7.8907e-01,-8.2077e-01,-3.8747e-01,-4.7485e-01,-8.4246e-01,-6.5592e-01;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
	{
		const double p = 20.;
		const double power = p/(p-1.);
		LpDualityMapping J_p(p);
		Eigen::MatrixXd expected(2,10);
		expected << 9.8144e-10,-9.8730e-01,1.8154e-13,9.3566e-01,1.1623e-19,-5.7242e-06,-9.5207e-01,7.7562e-01,-1.4377e-18,-6.1100e-02,
				9.7424e-01,1.4816e-03,9.3184e-01,-1.6496e-02,-9.5135e-01,-9.6554e-01,-4.8824e-02,-2.3461e-01,-9.6456e-01,-9.2434e-01;
		for (size_t i=0; i<10; ++i) {
			const Eigen::VectorXd compare = J_p(X.col(i), power);
//			std::cout << "# " << i << ": Expecting " << expected.col(i).transpose()
//					<< " and got " << compare.transpose() << ".\n";
			CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
		}
	}
}

void LpDualityMappingUnitTest::setTolerance()
{
	LpDualityMapping J_infty(0);
	CPPUNIT_ASSERT_EQUAL(BASSOTOLERANCE, J_infty.tolerance);
	const double value = 1e-1;
	J_infty.setTolerance(value);
	CPPUNIT_ASSERT_EQUAL(value, J_infty.tolerance);
}
