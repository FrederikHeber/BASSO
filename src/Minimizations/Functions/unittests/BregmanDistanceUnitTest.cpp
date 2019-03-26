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
 * BregmanDistanceUnitTest.cpp
 *
 *  Created on: Apr 25, 2014
 *      Author: heber
 */

#include "BregmanDistanceUnitTest.hpp"

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <Eigen/Dense>
#include <limits>

#include "Log/Logging.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Functions/BregmanDistance.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( BregmanDistanceUnitTest );

using namespace boost::assign;

void BregmanDistanceUnitTest::setUp()
{
	// BregmanDistance uses logging
	boost::log::core::get()->set_filter
			(
					boost::log::trivial::severity >= boost::log::trivial::info
			);
	startLogging();
}


void BregmanDistanceUnitTest::tearDown()
{
}

/** We generate test vector as follows via octave:
 *
 * function testvectors(p,q)
 * 	x=ones(2,10)-2.*rand(2,10)
 * 	y=ones(2,10)-2.*rand(2,10)
 * 	for i=1:10
 * 		fval=BregmanDistance(x(:,i),y(:,i),p,q,1e-6)
 * 	endfor
 * endfunction
 *
 * with
 *
 * function fval=BregmanDistance(x,y,p,power,Tol)
 * 	q=power/(power-1);
 * 	fval=1/q*norm(x,p)^power;
 * 	fval+=1/power*norm(y,p)^power;
 * 	fval-=DualityMapping(x,p,power,Tol)'*y;
 * endfunction
 */

void BregmanDistanceUnitTest::oneoneNorm()
{
	const double p = 1.1;
	Eigen::MatrixXd X(2,10);
	Eigen::MatrixXd Y(2,10);
	X << -0.921969,-0.023463,0.879205,-0.085334,0.075672,0.712906,-0.643552,-0.996276,0.676741,-0.937033,
			-0.852697,0.657800,0.468902,0.842803,-0.130822,0.896841,-0.616836,-0.570879,-0.721382,-0.864484;
	Y << 0.540338,-0.535169,0.942398,0.640759,-0.784199,0.709511,0.500785,0.227726,0.134048,-0.059655,
			-0.315758,-0.139483,-0.522375,0.069016,-0.743665,0.111088,0.451357,0.239402,0.841339,0.548371;
	{
		const double power = 1.01;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 1.0286,  0.39996,  0.94789,  1.1274,  1.4263,  0.031599,  1.7916,  0.88372,  1.6141,  1.0585;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 1.1;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 1.1023,  0.38601,  0.96904,  1.1197,  1.3661,  0.051830,  1.8225,  0.96079,  1.6595,  1.1536;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 2.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 2.0704,  0.27066,  1.2085,  1.0431,  1.0478,  0.30330,  2.1581,  1.8317,  2.1812,  2.3781;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
}

void BregmanDistanceUnitTest::onefiveNorm()
{
	const double p = 1.5;
	Eigen::MatrixXd X(2,10);
	Eigen::MatrixXd Y(2,10);
	X << -0.921969,-0.023463,0.879205,-0.085334,0.075672,0.712906,-0.643552,-0.996276,0.676741,-0.937033,
			-0.852697,0.657800,0.468902,0.842803,-0.130822,0.896841,-0.616836,-0.570879,-0.721382,-0.864484;
	Y << 0.540338,-0.535169,0.942398,0.640759,-0.784199,0.709511,0.500785,0.227726,0.134048,-0.059655,
			-0.315758,-0.139483,-0.522375,0.069016,-0.743665,0.111088,0.451357,0.239402,0.841339,0.548371;
	{
		const double power = 1.1;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 0.93625,  0.59528,  0.69043,  0.78009,  1.0278,  0.13311,  1.5156,  0.79453,  1.4681,  1.0083;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 1.5;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 1.1679,  0.50592,  0.71748,  0.74407,  0.85996,  0.20259,  1.5282,  1.0340,  1.5413,  1.3028;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 2;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 1.5014,  0.41277,  0.75285,  0.70039,  0.72770,  0.29689,  1.5428,  1.3305,  1.6369,  1.7214;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
}

void BregmanDistanceUnitTest::twoNorm()
{
	const double p = 2.;
	Eigen::MatrixXd X(2,10);
	Eigen::MatrixXd Y(2,10);
	X << -0.921969,-0.023463,0.879205,-0.085334,0.075672,0.712906,-0.643552,-0.996276,0.676741,-0.937033,
			-0.852697,0.657800,0.468902,0.842803,-0.130822,0.896841,-0.616836,-0.570879,-0.721382,-0.864484;
	Y << 0.540338,-0.535169,0.942398,0.640759,-0.784199,0.709511,0.500785,0.227726,0.134048,-0.059655,
			-0.315758,-0.139483,-0.522375,0.069016,-0.743665,0.111088,0.451357,0.239402,0.841339,0.548371;
	{
		const double power = 1.5;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 1.0035,  0.54981,  0.49253,  0.60101,  0.67100,  0.24884,  1.2858,  0.87602,  1.3712,  1.1233;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 2.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 1.2133,  0.44875,  0.49331,  0.56298,  0.55748,  0.30871,  1.2253,  1.0774,  1.3683,  1.3830;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 4.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 2.2649,  0.19848,  0.49686,  0.42683,  0.34058,  0.56393,  1.0026,  1.7861,  1.3547,  2.6842;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
}

void BregmanDistanceUnitTest::threeNorm()
{
	const double p = 3.;
	Eigen::MatrixXd X(2,10);
	Eigen::MatrixXd Y(2,10);
	X << -0.921969,-0.023463,0.879205,-0.085334,0.075672,0.712906,-0.643552,-0.996276,0.676741,-0.937033,
			-0.852697,0.657800,0.468902,0.842803,-0.130822,0.896841,-0.616836,-0.570879,-0.721382,-0.864484;
	Y << 0.540338,-0.535169,0.942398,0.640759,-0.784199,0.709511,0.500785,0.227726,0.134048,-0.059655,
			-0.315758,-0.139483,-0.522375,0.069016,-0.743665,0.111088,0.451357,0.239402,0.841339,0.548371;
	{
		const double power = 1.5;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 0.87883,  0.55372,  0.26112,  0.54290,  0.48807,  0.31396,  1.0825,  0.74834,  1.2462,  0.96980;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 3.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 1.2286,  0.30182,  0.23470,  0.44296,  0.29139,  0.39202,  0.78581,  1.0958,  1.0326,  1.3917;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 10.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 3.2942,  0.017059,  0.14461,  0.15105,  0.068712,  0.63698,  0.16631,  1.9832,  0.42931,  4.1186;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
}

void BregmanDistanceUnitTest::sixNorm()
{
	const double p = 6.;
	Eigen::MatrixXd X(2,10);
	Eigen::MatrixXd Y(2,10);
	X << -0.921969,-0.023463,0.879205,-0.085334,0.075672,0.712906,-0.643552,-0.996276,0.676741,-0.937033,
			-0.852697,0.657800,0.468902,0.842803,-0.130822,0.896841,-0.616836,-0.570879,-0.721382,-0.864484;
	Y << 0.540338,-0.535169,0.942398,0.640759,-0.784199,0.709511,0.500785,0.227726,0.134048,-0.059655,
			-0.315758,-0.139483,-0.522375,0.069016,-0.743665,0.111088,0.451357,0.239402,0.841339,0.548371;
	{
		const double power = 2.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 0.86548,  0.45132,  0.041182,  0.50228,  0.28914,  0.42560,  0.77550,  0.77262,  1.0428,  0.87435;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 6.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 1.0541,  0.088607,  0.030654,  0.28085,  0.066931,  0.36918,  0.20473,  1.0818,  0.40193,  1.1381;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 10.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_p = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_p,
				J_p.getPower());
		Eigen::VectorXd expected(10);
		expected << 1.1155,  0.017061,  0.022774,  0.14910,  0.021871,  0.29738,  0.052789,  1.1588,  0.15511,  1.2815;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
}


void BregmanDistanceUnitTest::inftyNorm()
{
	const double p = std::numeric_limits<double>::infinity();
	Eigen::MatrixXd X(2,10);
	Eigen::MatrixXd Y(2,10);
	X << -0.921969,-0.023463,0.879205,-0.085334,0.075672,0.712906,-0.643552,-0.996276,0.676741,-0.937033,
			-0.852697,0.657800,0.468902,0.842803,-0.130822,0.896841,-0.616836,-0.570879,-0.721382,-0.864484;
	Y << 0.540338,-0.535169,0.942398,0.640759,-0.784199,0.709511,0.500785,0.227726,0.134048,-0.059655,
			-0.315758,-0.139483,-0.522375,0.069016,-0.743665,0.111088,0.451357,0.239402,0.841339,0.548371;
	{
		const double power = 1.1;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_infty = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_infty,
				J_infty.getPower());
		Eigen::VectorXd expected(10);
		expected << 1.0810,  0.64814,  2.1952e-04,  0.56462,  0.098692,  0.59401,  0.96001,  0.50682,  1.6295,  0.49481;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 2.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_infty = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_infty,
				J_infty.getPower());
		Eigen::VectorXd expected(10);
		expected << 1.0692,  0.45130,  0.0019967,  0.50228,  0.21875,  0.55424,  0.65475,  0.75182,  1.2210,  0.53347;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
	{
		const double power = 10.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(power);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(X.innerSize(), "lp", args);
		const Mapping &J_infty = *SpaceX->getSpace()->getDualityMapping();
		const BregmanDistance d_p(
				*SpaceX->getSpace()->getNorm(),
				J_infty,
				J_infty.getPower());
		Eigen::VectorXd expected(10);
		expected << 0.65969,  0.017060,  0.0078148,  0.14910,  0.0087956,  0.26450,  0.020548,  1.0873,  0.096628,  0.43669;
		for (size_t i=0; i<10; ++i) {
			SpaceElement_ptr_t x = ElementCreator::create(SpaceX,X.col(i));
			SpaceElement_ptr_t y = ElementCreator::create(SpaceX,Y.col(i));
//			std::cout << "# " << i << ": Expecting " << expected(i) << " and got " << d_p(x,y) << ".\n";
			CPPUNIT_ASSERT( fabs( (expected(i) - d_p(x,y) ) ) < 1e-4);
		}
	}
}
