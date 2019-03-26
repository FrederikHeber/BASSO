/*
 * point_tUnitTest.cpp
 *
 *  Created on: Jul 07, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/assign.hpp>
#include <ComputerTomography/DiscretizedRadon/unittests/point_tUnitTest.hpp>

#include <vector>

#include "Log/Logging.hpp"

#include "ComputerTomography/DiscretizedRadon/point_t.hpp"

using namespace boost::assign;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( point_tUnitTest );

void point_tUnitTest::setUp()
{
	// DiscretizedRadonMatrix uses logging
	boost::log::core::get()->set_filter
			(
					boost::log::trivial::severity >= boost::log::trivial::info
			);
	startLogging();
}


void point_tUnitTest::tearDown()
{
}

static point_t constructPoint(const double _x, const double _y)
{
	point_t point;
	point[0]=_x;
	point[1]=_y;
	return point;
}

void point_tUnitTest::getSquaredDistanceTest()
{
	typedef std::vector<point_t> points_t;
	points_t points;
	points +=
			constructPoint(0.,0.),
			constructPoint(0.1,0.15),
			constructPoint(1.,0.),
			constructPoint(1.2,1.4);

	std::vector<double> distances;
	distances +=
			0.0325,
			1.,
			3.4,
			0.8325,
			2.7725,
			2.;

	std::vector<double>::const_iterator distanceiter = distances.begin();
	for (points_t::const_iterator firstiter = points.begin();
			firstiter != points.end(); ++firstiter)
		for (points_t::const_iterator seconditer = firstiter;
				seconditer != points.end(); ++seconditer) {
			if (firstiter == seconditer)
				continue;
			const double distance = ((*firstiter) - (*seconditer)).squaredNorm();
//			std::cout << "True distance is " << distance << " for points "
//					<< *firstiter << " and " << *seconditer
//					<< ", expecting " << *distanceiter << std::endl;
			CPPUNIT_ASSERT( fabs(distance - *distanceiter) < BASSOTOLERANCE );
			++distanceiter;
		}
	CPPUNIT_ASSERT( distanceiter == distances.end() );
}

void point_tUnitTest::maxNormTest()
{
	typedef std::vector<point_t> points_t;
	points_t points;
	points +=
			constructPoint(0.,0.),
			constructPoint(0.1,-0.15),
			constructPoint(-1.,0.),
			constructPoint(1.2,-1.4);

	std::vector<double> norms;
	norms +=
			0.,
			0.15,
			1.,
			1.4;

	std::vector<double>::const_iterator normiter = norms.begin();
	for (points_t::const_iterator iter = points.begin();
			iter != points.end(); ++iter) {
		const double maxnorm = (*iter).lpNorm<Eigen::Infinity>();
		CPPUNIT_ASSERT( fabs(maxnorm - *normiter) < BASSOTOLERANCE );
		++normiter;
	}
}
