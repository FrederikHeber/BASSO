/*
 * DiscretizedRadonMatrixUnitTest.cpp
 *
 *  Created on: Jul 07, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include "ComputerTomographyBase/DiscretizedRadon/unittests/DiscretizedRadonMatrixUnitTest.hpp"

#include <Eigen/Dense>

#include "Log/Logging.hpp"

#include <ComputerTomographyBase/DiscretizedRadon/BackprojectionMatrix.hpp>
#include "ComputerTomographyBase/DiscretizedRadon/DiscretizedRadonMatrix.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DiscretizedRadonMatrixUnitTest );

void DiscretizedRadonMatrixUnitTest::setUp()
{
	// DiscretizedRadonMatrix uses logging
	logging::core::get()->set_filter
			(
					logging::trivial::severity >= logging::trivial::debug
			);
	startLogging();
}


void DiscretizedRadonMatrixUnitTest::tearDown()
{
}

const unsigned int num_angles = 4;
const unsigned int num_offsets = 5;
const unsigned int num_pixels = 5;

void DiscretizedRadonMatrixUnitTest::adjointTest()
{
	Eigen::VectorXd f(num_pixels*num_pixels);
	f.setRandom();
	Eigen::VectorXd g(num_angles*num_offsets);
	g.setRandom();
	DiscretizedRadonMatrix RadonMatrix(num_pixels, num_pixels, num_angles, num_offsets);
	BackprojectionMatrix backprojection(num_pixels, num_pixels, num_angles, num_offsets);
	// testing for adjoint correctness
	{
		const Eigen::MatrixXd& R = RadonMatrix.getMatrix();
		const Eigen::MatrixXd Rs = RadonMatrix.getMatrix().transpose();
		const Eigen::VectorXd Rf = R*f;
		const Eigen::VectorXd Rsg = Rs*g;

		std::cout << "Got f=" << f.transpose() << std::endl;
		std::cout << "Got Rf=" << Rf.transpose() << std::endl;
		std::cout << "Got g=" << g.transpose() << std::endl;
		std::cout << "Got Rsg=" << Rsg.transpose() << std::endl;
		std::cout << "Then, <Rf,g> = " << g.dot(Rf)
				<< " and <f,Rsg> = " << f.dot(Rsg) << std::endl;

		CPPUNIT_ASSERT( fabs(g.dot(Rf) - f.dot(Rsg)) < BASSOTOLERANCE );
	}

	{
		const Eigen::MatrixXd& Rs = backprojection.getMatrix();
		const Eigen::MatrixXd R = backprojection.getMatrix().transpose();
		const Eigen::VectorXd Rf = R*f;
		const Eigen::VectorXd Rsg = Rs*g;

		std::cout << "Got f=" << f.transpose() << std::endl;
		std::cout << "Got Rf=" << Rf.transpose() << std::endl;
		std::cout << "Got g=" << g.transpose() << std::endl;
		std::cout << "Got Rsg=" << Rsg.transpose() << std::endl;
		std::cout << "Then, <Rf,g> = " << g.dot(Rf)
				<< " and <f,Rsg> = " << f.dot(Rsg) << std::endl;

		CPPUNIT_ASSERT( fabs(g.dot(Rf) - f.dot(Rsg)) < BASSOTOLERANCE );
	}

	{

		const Eigen::MatrixXd& R = RadonMatrix.getMatrix();
		const Eigen::MatrixXd Rs = backprojection.getMatrix();
		const Eigen::VectorXd Rf = R*f;
		const Eigen::VectorXd Rsg = Rs*g;

		std::cout << "Got f=" << f.transpose() << std::endl;
		std::cout << "Got Rf=" << Rf.transpose() << std::endl;
		std::cout << "Got g=" << g.transpose() << std::endl;
		std::cout << "Got Rsg=" << Rsg.transpose() << std::endl;
		std::cout << "Then, <Rf,g> = " << g.dot(Rf)
				<< " and <f,Rsg> = " << f.dot(Rsg) << std::endl;

		std::cout << "R=" << R << std::endl;
		std::cout << "Rs'=" << Rs.transpose() << std::endl;
		std::cout << R-Rs.transpose() << std::endl;

		// The following should work only for true adjoint
//		CPPUNIT_ASSERT( fabs(g.dot(Rf) - f.dot(Rsg)) < BASSOTOLERANCE );
	}
}
