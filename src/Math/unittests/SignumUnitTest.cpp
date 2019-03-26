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
 * SignumUnitTest.cpp
 *
 *  Created on: Jul 14, 2014
 *      Author: heber
 */

#include "SignumUnitTest.hpp"

#include <Eigen/Dense>

#include "Math/Helpers.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SignumUnitTest );


void SignumUnitTest::setUp()
{
}


void SignumUnitTest::tearDown()
{
}

/** We generate test vectors as follows via octave:
 *
 * x=ones(2,10)-2.*rand(2,10)
 * sign(x)
 *
 */


void SignumUnitTest::twoTest()
{
	Eigen::MatrixXd X(2,10);
	X << 0.762452,0.749991,-0.823097,-0.649048,0.232539,0.824522,0.496256,0.659914,0.700726,-0.519157,
			0.388856,0.012850,-0.914217,0.551940,0.827947,-0.245284,-0.345104,-0.143636,-0.438064,-0.448346;
	Eigen::MatrixXd expected(2,10);
	expected << 1,1,-1,-1,1,1,1,1,1,-1,
			1,1,-1,1,1,-1,-1,-1,-1,-1;
	for (size_t i=0; i<10; ++i) {
		const Eigen::VectorXd compare = Helpers::signum(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
		CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
	}
}

void SignumUnitTest::fiveTest()
{
	Eigen::MatrixXd X(5,10);
	X << -0.0510168,-0.6125603,0.8096526,-0.4378988,-0.1701139,-0.8221936,-0.4176625,-0.4524001,-0.3210781,-0.4809454,
		0.8443800,0.3118541,-0.2123250,0.9929333,-0.6989529,-0.9901448,-0.6172886,-0.0990839,0.1191296,-0.5384438,
		0.5773228,-0.0087454,-0.8525228,0.1900142,0.2773021,-0.1277766,0.6077737,-0.3347506,0.1377782,-0.0893297,
		0.0088735,0.7995322,0.3846394,-0.6339381,0.8413644,0.7822711,-0.4677034,-0.5352645,-0.3219619,0.0626028,
		-0.9483627,-0.0231827,-0.4805136,-0.3156587,0.9363985,0.0834479,-0.6988783,0.2773671,-0.2777026,0.5446490;
	Eigen::MatrixXd expected(5,10);
	expected << -1,-1,1,-1,-1,-1,-1,-1,-1,-1,
		1,1,-1,1,-1,-1,-1,-1,1,-1,
		1,-1,-1,1,1,-1,1,-1,1,-1,
		1,1,1,-1,1,1,-1,-1,-1,1,
		-1,-1,-1,-1,1,1,-1,1,-1,1;
	for (size_t i=0; i<10; ++i) {
		const Eigen::VectorXd compare = Helpers::signum(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
		CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
	}
}

void SignumUnitTest::tenTest()
{
	Eigen::MatrixXd X(10,10);
	X << -0.2748037,0.1297174,0.0971120,-0.7302618,0.5159538,0.6696792,0.6163268,0.1666823,-0.9158217,0.4139321,
		-0.0060395,0.9082336,0.8451014,-0.1954476,-0.3341592,0.3618202,-0.4846105,-0.5535242,0.2616558,0.9519939,
		-0.0087249,-0.3982453,0.0752541,0.0685472,-0.1576404,-0.4396471,-0.3602312,-0.2463073,-0.5602789,0.5816376,
		-0.1160532,0.7387554,0.2043750,-0.5911679,-0.4468587,0.4944709,-0.6856535,-0.2229509,-0.6077783,0.0579071,
		0.0441697,-0.2136574,0.9970342,0.5330651,-0.0956938,-0.1833911,0.7226938,-0.0726103,-0.7399354,0.4597379,
		0.6422679,0.6519052,-0.1210018,-0.9680343,0.7073513,0.9660143,0.5244329,0.8021393,0.3526311,-0.8798049,
		0.3437433,-0.1041717,0.3336555,-0.5814434,-0.7391993,-0.7770869,0.9149192,0.8991434,0.2635078,-0.0360407,
		-0.4899371,0.2894166,-0.0367166,-0.6400291,-0.6298687,0.9090958,-0.1192517,0.3251074,-0.6573007,0.9469398,
		-0.2828731,0.1997232,-0.9090181,-0.5099530,-0.2862054,-0.7003964,-0.1950223,-0.4725289,0.3938916,0.1149939,
		-0.5577270,0.6559519,-0.8756973,-0.8727486,0.1669378,0.8277767,-0.4541044,0.1214144,-0.1868776,0.0570731;
	Eigen::MatrixXd expected(10,10);
	expected << -1,1,1,-1,1,1,1,1,-1,1,
		-1,1,1,-1,-1,1,-1,-1,1,1,
		-1,-1,1,1,-1,-1,-1,-1,-1,1,
		-1,1,1,-1,-1,1,-1,-1,-1,1,
		1,-1,1,1,-1,-1,1,-1,-1,1,
		1,1,-1,-1,1,1,1,1,1,-1,
		1,-1,1,-1,-1,-1,1,1,1,-1,
		-1,1,-1,-1,-1,1,-1,1,-1,1,
		-1,1,-1,-1,-1,-1,-1,-1,1,1,
		-1,1,-1,-1,1,1,-1,1,-1,1;
	for (size_t i=0; i<10; ++i) {
		const Eigen::VectorXd compare = Helpers::signum(X.col(i));
//			std::cout << "# " << i << ": Expecting " << expected.row(i)
//					<< " and got " << compare << ".\n";
		CPPUNIT_ASSERT( expected.col(i).isApprox(compare, 1e-4)  );
	}
}
