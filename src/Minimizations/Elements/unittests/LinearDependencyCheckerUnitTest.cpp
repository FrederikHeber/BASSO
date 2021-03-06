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
 * LinearDependencyCheckerUnitTest.cpp
 *
 *  Created on: Feb 10, 2015
 *      Author: heber
 */

#include <boost/assign.hpp>
#include <Minimizations/Elements/unittests/LinearDependencyCheckerUnitTest.hpp>
#include <iostream>

#include "Log/Logging.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

using namespace boost::assign;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( LinearDependencyCheckerUnitTest );

void LinearDependencyCheckerUnitTest::setUp()
{
	// LinearDependencyChecker uses logging
	boost::log::core::get()->set_filter
			(
					boost::log::trivial::severity >= boost::log::trivial::info
			);
	startLogging();
}


void LinearDependencyCheckerUnitTest::tearDown()
{
}

void LinearDependencyCheckerUnitTest::noVector()
{
	std::cout << "The following warning is desired and does not imply a failure of the test.\n";
	LinearDependencyChecker::vectors_t vectors(0);
	CPPUNIT_ASSERT( lindepcheck(vectors) );
}

void LinearDependencyCheckerUnitTest::singleVector()
{
	const unsigned int dim = 5;
	const double p=2.;
	NormedSpaceFactory::args_t args;
	args += boost::any(p), boost::any(p);
	const NormedSpace_ptr_t SpaceX =
			NormedSpaceFactory::create(dim, "lp", args);

	// single zero vector
	{
		SpaceElement_ptr_t x = SpaceX->createElement();
		x->setZero();
		LinearDependencyChecker::vectors_t vectors(1, x);
		CPPUNIT_ASSERT( lindepcheck(vectors) );
	}

	// single non zero vector
	{
		Eigen::VectorXd xaxis(dim);
		xaxis.setZero();
		xaxis[0] = 1;
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX, xaxis);
		LinearDependencyChecker::vectors_t vectors(1, x);
		CPPUNIT_ASSERT( !lindepcheck(vectors) );
	}
}


void LinearDependencyCheckerUnitTest::twoVectors()
{
	const unsigned int dim = 5;
	const double p=2.;
	NormedSpaceFactory::args_t args;
	args += boost::any(p), boost::any(p);
	const NormedSpace_ptr_t SpaceX =
			NormedSpaceFactory::create(dim, "lp", args);

	// two zero vectors
	{
		SpaceElement_ptr_t x = SpaceX->createElement();
		x->setZero();
		SpaceElement_ptr_t y = SpaceX->createElement();
		y->setZero();
		LinearDependencyChecker::vectors_t vectors;
		vectors += x,y;
		CPPUNIT_ASSERT( lindepcheck(vectors) );
	}

	// single non zero vector and zero vector
	{
		Eigen::VectorXd xaxis(dim);
		xaxis.setZero();
		xaxis[0] = 1;
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX, xaxis);
		SpaceElement_ptr_t y = SpaceX->createElement();
		y->setZero();
		LinearDependencyChecker::vectors_t vectors;
		vectors += x,y;
		CPPUNIT_ASSERT( lindepcheck(vectors) );
	}

	// two non-zero but dependent vectors
	{
		Eigen::VectorXd xaxis(dim);
		xaxis.setZero();
		xaxis[0] = 1;
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX, xaxis);
		SpaceElement_ptr_t y = ElementCreator::create(SpaceX, xaxis);
		LinearDependencyChecker::vectors_t vectors;
		vectors += x,y;
		CPPUNIT_ASSERT( lindepcheck(vectors) );
	}

	// two non-zero independent orthogonalvectors
	{
		Eigen::VectorXd xaxis(dim);
		xaxis.setZero();
		xaxis[0] = 1;
		Eigen::VectorXd yaxis(dim);
		yaxis.setZero();
		yaxis[1] = 1;
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX, xaxis);
		SpaceElement_ptr_t y = ElementCreator::create(SpaceX, yaxis);
		LinearDependencyChecker::vectors_t vectors;
		vectors += x,y;
		CPPUNIT_ASSERT( !lindepcheck(vectors) );
	}

	// two non-zero complex independent orthogonal vectors
	{
		Eigen::VectorXd xaxis(dim);
		xaxis.setZero();
		xaxis[0] = 1./sqrt(2.);
		xaxis[1] = 1./sqrt(2.);
		Eigen::VectorXd yaxis(dim);
		yaxis.setZero();
		yaxis[0] = -1./sqrt(2.);
		yaxis[1] = 1./sqrt(2.);
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX, xaxis);
		SpaceElement_ptr_t y = ElementCreator::create(SpaceX, yaxis);
		LinearDependencyChecker::vectors_t vectors;
		vectors += x,y;
		CPPUNIT_ASSERT( !lindepcheck(vectors) );
	}

	// two non-zero complex independent non-orthogonal vectors
	{
		Eigen::VectorXd xaxis(dim);
		xaxis.setZero();
		xaxis[0] = 0.234;
		xaxis[1] = 0.5345;
		xaxis[2] = 0.54555;
		Eigen::VectorXd yaxis(dim);
		yaxis.setZero();
		yaxis[2] = 2.34782;
		yaxis[3] = 0.74563;
		yaxis[4] = 3.8564;
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX, xaxis);
		SpaceElement_ptr_t y = ElementCreator::create(SpaceX, yaxis);
		LinearDependencyChecker::vectors_t vectors;
		vectors += x,y;
		CPPUNIT_ASSERT( !lindepcheck(vectors) );
	}
}

void LinearDependencyCheckerUnitTest::threeVectors()
{
	// three vectors of dim two
	{
		const unsigned int dim = 2;
		const double p=2.;
		NormedSpaceFactory::args_t args;
		args += boost::any(p), boost::any(p);
		const NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::create(dim, "lp", args);

		Eigen::VectorXd xaxis(dim);
		xaxis.setZero();
		xaxis[0] = 1.;
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX, xaxis);
		Eigen::VectorXd yaxis(dim);
		yaxis.setZero();
		yaxis[1] = 1.;
		SpaceElement_ptr_t y = ElementCreator::create(SpaceX, yaxis);
		Eigen::VectorXd xyaxis(dim);
		xyaxis.setZero();
		xyaxis[0] = -1./sqrt(2.);
		xyaxis[1] = 1./sqrt(2.);
		SpaceElement_ptr_t xy = ElementCreator::create(SpaceX, xyaxis);
		LinearDependencyChecker::vectors_t vectors;
		vectors += x,y,xy;
		CPPUNIT_ASSERT( lindepcheck(vectors) );
	}

	const unsigned int dim = 5;
	const double p=2.;
	NormedSpaceFactory::args_t args;
	args += boost::any(p), boost::any(p);
	const NormedSpace_ptr_t SpaceX =
			NormedSpaceFactory::create(dim, "lp", args);

	// three non-zero complex independent non-orthogonal vectors
	{
		Eigen::VectorXd xaxis(dim);
		xaxis.setZero();
		xaxis[0] = 0.234;
		xaxis[1] = 0.5345;
		xaxis[2] = 0.54555;
		Eigen::VectorXd yaxis(dim);
		yaxis.setZero();
		yaxis[2] = 2.34782;
		yaxis[3] = 0.74563;
		yaxis[4] = 3.8564;
		Eigen::VectorXd zaxis(dim);
		zaxis.setZero();
		zaxis[1] = 1.857604782;
		zaxis[2] = 1.122563;
		zaxis[3] = 0.3464;
		SpaceElement_ptr_t x = ElementCreator::create(SpaceX, xaxis);
		SpaceElement_ptr_t y = ElementCreator::create(SpaceX, yaxis);
		SpaceElement_ptr_t z = ElementCreator::create(SpaceX, zaxis);
		LinearDependencyChecker::vectors_t vectors;
		vectors += x,y,z;
		CPPUNIT_ASSERT( !lindepcheck(vectors) );
	}
}
