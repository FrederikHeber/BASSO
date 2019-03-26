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
		expected <<
				0.5478745260841271,
				1.541107636450696,
				0.9424331396506981,
				0.6045384472552063,
				0.2391841275896053,
				1.000046882477412,
				0.9204524054224361,
				1.629028269105949,
				0.5183939146105486,
				0.7288677395775812;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				0.5278650127819026,
				1.304548400199131,
				0.8491539867770614,
				0.5270252446811835,
				0.2045598746557526,
				0.9360718117533751,
				0.8902991206445896,
				1.383326740562306,
				0.4438595082991083,
				0.6424765422467946;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				0.5208780383182567,
				0.9420680723484283,
				0.780144945860734,
				0.4471827896676158,
				0.1594871240480309,
				0.9037437271937642,
				0.8804180330299162,
				1.026914689852702,
				0.3483657964429822,
				0.5632508552077028;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				0.5223076852249064,
				1.164233679492652,
				0.8068206901939241,
				0.4861045919408702,
				0.1848128563601569,
				0.9129487046411753,
				0.8823223801111474,
				1.239313726626555,
				0.401540038941324,
				0.5989906666551659;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				1.379204618735429,
				1.762807446234149,
				1.449916341854767,
				0.4631167945397467,
				1.530356674825692,
				1.043928877980617,
				0.8323476379124169,
				1.674116849265731,
				0.8195304746886963,
				1.183676001026054;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				1.897584096427657,
				2.136648110246748,
				1.546168795024981,
				1.568028070746475,
				1.922045574180019,
				1.426097313171861,
				1.732441689374854,
				1.363296056247872,
				2.132786287926899,
				1.90067685274299;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				0.5209487539672338,
				1.040770548963172,
				0.7848319642513533,
				0.4580712862990062,
				0.1689258944084371,
				0.9046230945798833,
				0.8805016320324537,
				1.116078030524289,
				0.3678248372808625,
				0.5718058616779829;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				1.156925308143252,
				1.380591110982202,
				1.190750108539649,
				0.3862377363534127,
				1.240214316762886,
				0.872336461239226,
				0.685505044725172,
				1.316485337885209,
				0.7025531281504473,
				1.007544317069905;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				1.408278528885601,
				1.565617155958245,
				1.16091195702887,
				1.26525692537547,
				1.448039621498664,
				1.110886445079814,
				1.275803356441185,
				1.05287465727844,
				1.563611976113804,
				1.386910628509519;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				0.5208819348767975,
				0.9856564499812991,
				0.7810173538555848,
				0.4502902356827378,
				0.1630029514643971,
				0.9038357554111797,
				0.8804221281378948,
				1.064151591495952,
				0.3554898803089913,
				0.5653675324332494;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				1.073015242795003,
				1.228870112521734,
				1.086030516798965,
				0.3574294629077092,
				1.127952682305419,
				0.8053536407350824,
				0.6323609050903276,
				1.179960009295014,
				0.6608732054435563,
				0.949295881450161;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				1.231436302119236,
				1.354087159309278,
				1.036482455736694,
				1.150415961989,
				1.276534889888392,
				0.9984955580940865,
				1.111543552178988,
				0.9393945992257847,
				1.354710889299011,
				1.201354126192843;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				0.5208780000000176,
				0.9008391471410394,
				0.7800700309472804,
				0.4465128080603014,
				0.157887505566148,
				0.9037410000423686,
				0.880418000000008,
				1.000509317084285,
				0.3453051005190296,
				0.5629227529945779;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				0.9793351851895178,
				1.0089615415637,
				0.9277769799143303,
				0.3214918161871571,
				0.9736647972208938,
				0.7144815677208166,
				0.5849836575175007,
				1.009564053265138,
				0.6211794508429107,
				0.9011251293104742;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
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
		expected <<
				1.013255386054021,
				1.053374394701754,
				0.9262779376510162,
				1.00311790343948,
				1.043004453896559,
				0.8568774326775785,
				0.9026249276205977,
				0.7900631268082297,
				1.072090314459995,
				0.9861739728219595;

		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX, X.col(i));
			const double compare = norm(x);
//			std::cout << std::setprecision(16) << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
			CPPUNIT_ASSERT( fabs( expected(i,0) - compare) < BASSOTOLERANCE  );
		}
	}
}
