/*
 * BackprojectionMatrixUnitTest.cpp
 *
 *  Created on: Jul 07, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include "ComputerTomographyBase/DiscretizedRadon/unittests/BackprojectionMatrixUnitTest.hpp"

#include <Eigen/Dense>

#include "Log/Logging.hpp"

#include <ComputerTomographyBase/DiscretizedRadon/BackprojectionMatrix.hpp>
#include "ComputerTomographyBase/DiscretizedRadon/BackprojectionMatrix.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( BackprojectionMatrixUnitTest );

void BackprojectionMatrixUnitTest::setUp()
{
	// BackprojectionMatrix uses logging
	boost::log::core::get()->set_filter
			(
					boost::log::trivial::severity >= boost::log::trivial::debug
			);
	startLogging();
}


void BackprojectionMatrixUnitTest::tearDown()
{
}

void BackprojectionMatrixUnitTest::singlePixel_few_angles_few_offsetsTest()
{
	// constants
	const int num_pixel_x = 1;
	const int num_pixel_y = 1;
	const int num_pixels = num_pixel_x*num_pixel_y;

	const int num_offsets = 1;
	const int num_angles = 2;
	const int num_measurements = num_angles*num_offsets;

	// construct matrix
	BackprojectionMatrix Backprojection(
			num_pixel_x,
			num_pixel_y,
			num_angles,
			num_offsets);
	const Eigen::MatrixXd &A = Backprojection.getMatrix();

	// single pixel with white value
	{
		// create expected value
		Eigen::VectorXd x(num_measurements);
		x.setZero();
		x(0) = 1.;
		x(1) = 1.;

		Eigen::VectorXd expected_y(num_pixels);
		expected_y.setZero();
		expected_y(0) = (1.+1.)/2.*M_PI;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}
	// single pixel with grey value
	{
		// create expected value
		Eigen::VectorXd x(num_measurements);
		x.setZero();
		x(0) = .5;
		x(1) = .5;

		Eigen::VectorXd expected_y(num_pixels);
		expected_y.setZero();
		expected_y(0) = (.5+.5)/2.*M_PI;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}
}

void BackprojectionMatrixUnitTest::singlePixel_many_angles_few_offsetsTest()
{
	// constants
	const int num_pixel_x = 1;
	const int num_pixel_y = 1;
	const int num_pixels = num_pixel_x*num_pixel_y;

	const int num_offsets = 1;
	const int num_angles = 8;
	const int num_measurements = num_angles*num_offsets;

	// construct matrix
	BackprojectionMatrix Backprojection(
			num_pixel_x,
			num_pixel_y,
			num_angles,
			num_offsets);
	const Eigen::MatrixXd &A = Backprojection.getMatrix();

	// single pixel with white value one
	{
		// create expected value
		Eigen::VectorXd x(num_measurements);
		x.setZero();
		x(0) = 1.;
		x(1) = 1.;
		x(2) = 1.;
		x(3) = 1.;
		x(4) = 1.;
		x(5) = 1.;
		x(6) = 1.;
		x(7) = 1.;

		Eigen::VectorXd expected_y(num_pixels);
		expected_y.setZero();
		expected_y(0) = (1.+1.+1.+1.+1.+1.+1.+1.)/8.*M_PI;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y);
	}
	// single pixel with grey value one
	{
		// create expected value
		Eigen::VectorXd x(num_measurements);
		x.setZero();
		x(0) = .5;
		x(1) = .5;
		x(2) = .5;
		x(3) = .5;
		x(4) = .5;
		x(5) = .5;
		x(6) = .5;
		x(7) = .5;

		Eigen::VectorXd expected_y(num_pixels);
		expected_y.setZero();
		expected_y(0) = (.5+.5+.5+.5+.5+.5+.5+.5)/8.*M_PI;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y);
	}
}

void BackprojectionMatrixUnitTest::singlePixel_few_angles_many_offsetsTest()
{
	// constants
	const int num_pixel_x = 1;
	const int num_pixel_y = 1;
	const int num_pixels = num_pixel_x*num_pixel_y;

	const int num_offsets = 3;
	const int num_angles = 2;
	const int num_measurements = num_angles*num_offsets;

	// construct matrix
	BackprojectionMatrix Backprojection(
			num_pixel_x,
			num_pixel_y,
			num_angles,
			num_offsets);
	const Eigen::MatrixXd &A = Backprojection.getMatrix();

	// single pixel with white value one
	{
		// create expected value
		Eigen::VectorXd x(num_measurements);
		x.setZero();
		x(0) = 1.;
		x(1) = 1.;
		x(2) = 1.;
		x(3) = 1.;
		x(4) = 1.;
		x(5) = 1.;

		Eigen::VectorXd expected_y(num_pixels);
		expected_y.setZero();
		expected_y(0) = 1./1.*M_PI;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y);
	}
	// single pixel with grey value one
	{
		// create expected value
		Eigen::VectorXd x(num_measurements);
		x.setZero();
		x(0) = .5;
		x(1) = .5;
		x(2) = .5;
		x(3) = .5;
		x(4) = .5;
		x(5) = .5;

		Eigen::VectorXd expected_y(num_pixels);
		expected_y.setZero();
		expected_y(0) = .5/1.*M_PI;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y);
	}
}

void BackprojectionMatrixUnitTest::ninePixel_few_angles_few_offsetsTest()
{
	// constants
	const int num_pixel_x = 3;
	const int num_pixel_y = 3;
	const int num_offsets = 1;
	const int num_angles = 2;

	// derived constants
	const int num_pixels = num_pixel_x*num_pixel_y;
	const int num_measurements = num_angles*num_offsets;

	// construct matrix
	BackprojectionMatrix Backprojection(
			num_pixel_x,
			num_pixel_y,
			num_angles,
			num_offsets);
	const Eigen::MatrixXd &A = Backprojection.getMatrix();

	// nine pixel with single white pixel in center
	{
		Eigen::VectorXd x(num_measurements);
		x.setZero();
		x(0) = 1./3.;
		x(1) = 1./3.;

		// create expected value: intensity is constant from every
		// angle, hence image must be constant as well
		Eigen::VectorXd expected_y(num_pixels);
		expected_y.setZero();
		for (unsigned int i=0;i<num_pixel_x;++i)
			for (unsigned int j=0;j<num_pixel_x;++j)
				expected_y(i + (j * num_pixel_y)) =
						(1/3.+1/3.)/2.*M_PI;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y);
	}
}


void BackprojectionMatrixUnitTest::ninePixel_few_angles_many_offsetsTest()
{
	// constants
	const int num_pixel_x = 3;
	const int num_pixel_y = 3;
	const int num_offsets = 3;
	const int num_angles = 2;

	// derived constants
	const int num_pixels = num_pixel_x*num_pixel_y;
	const int num_measurements = num_angles*num_offsets;

	// construct matrix
	BackprojectionMatrix Backprojection(
			num_pixel_x,
			num_pixel_y,
			num_angles,
			num_offsets);
	const Eigen::MatrixXd &A = Backprojection.getMatrix();

	// nine pixel with single white pixel in center
	{
		Eigen::VectorXd x(num_measurements);
		x.setZero();
		x(0) = 0.;
		x(1) = 1./3.;
		x(2) = 0.;
		x(3) = 0.;
		x(4) = 1./3.;
		x(5) = 0.;

		// create expected value
		Eigen::VectorXd expected_y(num_pixels);
		expected_y.setZero();
		expected_y(0 + (0 * num_pixel_y)) =
				((2./3.*x(0)+1./3.*x(1))+(2./3.*x(3)+1./3.*x(4)))/2.*M_PI;
		expected_y(1 + (0 * num_pixel_y)) = 2./9.*M_PI;
		expected_y(2 + (0 * num_pixel_y)) = 1./9.*M_PI;
		expected_y(0 + (1 * num_pixel_y)) = 2./9.*M_PI;
		expected_y(1 + (1 * num_pixel_y)) = 3./9.*M_PI;
		expected_y(2 + (1 * num_pixel_y)) = 2./9.*M_PI;
		expected_y(0 + (2 * num_pixel_y)) = 1./9.*M_PI;
		expected_y(1 + (2 * num_pixel_y)) = 2./9.*M_PI;
		expected_y(2 + (2 * num_pixel_y)) = 1./9.*M_PI;

		// calculate y
		Eigen::VectorXd y = A * x;
//		CPPUNIT_ASSERT_EQUAL( expected_y, y);
		CPPUNIT_ASSERT( (expected_y - y).norm() < BASSOTOLERANCE);
	}

	// nine pixel with single white central horizontal line
	{
		Eigen::VectorXd x(num_measurements);
		x.setZero();
		x(0) = 0.;
		x(1) = 1.;
		x(2) = 0.;
		x(3) = 1./3.;
		x(4) = 1./3.;
		x(5) = 1./3.;

		// create expected value
		Eigen::VectorXd expected_y(num_pixels);
		expected_y.setZero();
		expected_y(0 + (0 * num_pixel_y)) =
				((2./3.*x(0)+1./3.*x(1))+(2./3.*x(3)+1./3.*x(4)))/2.*M_PI;
		expected_y(1 + (0 * num_pixel_y)) = 1./3.*M_PI;
		expected_y(2 + (0 * num_pixel_y)) = 1./3.*M_PI;
		expected_y(0 + (1 * num_pixel_y)) = 2./3.*M_PI;
		expected_y(1 + (1 * num_pixel_y)) = 2./3.*M_PI;
		expected_y(2 + (1 * num_pixel_y)) = 2./3.*M_PI;
		expected_y(0 + (2 * num_pixel_y)) = 1./3.*M_PI;
		expected_y(1 + (2 * num_pixel_y)) = 1./3.*M_PI;
		expected_y(2 + (2 * num_pixel_y)) = 1./3.*M_PI;

		// calculate y
		Eigen::VectorXd y = A * x;
//		CPPUNIT_ASSERT_EQUAL( expected_y, y);
		CPPUNIT_ASSERT( (expected_y - y).norm() < BASSOTOLERANCE);
	}

	// nine pixel with single white central vertical line
	{
		Eigen::VectorXd x(num_measurements);
		x.setZero();
		x(0) = 1./3.;
		x(1) = 1./3.;
		x(2) = 1./3.;
		x(3) = 0.;
		x(4) = 1.;
		x(5) = 0.;

		// create expected value
		Eigen::VectorXd expected_y(num_pixels);
		expected_y.setZero();
		expected_y(0 + (0 * num_pixel_y)) =
				((2./3.*x(0)+1./3.*x(1))+(2./3.*x(3)+1./3.*x(4)))/2.*M_PI;
		expected_y(1 + (0 * num_pixel_y)) = 2./3.*M_PI;
		expected_y(2 + (0 * num_pixel_y)) = 1./3.*M_PI;
		expected_y(0 + (1 * num_pixel_y)) = 1./3.*M_PI;
		expected_y(1 + (1 * num_pixel_y)) = 2./3.*M_PI;
		expected_y(2 + (1 * num_pixel_y)) = 1./3.*M_PI;
		expected_y(0 + (2 * num_pixel_y)) = 1./3.*M_PI;
		expected_y(1 + (2 * num_pixel_y)) = 2./3.*M_PI;
		expected_y(2 + (2 * num_pixel_y)) = 1./3.*M_PI;

		// calculate y
		Eigen::VectorXd y = A * x;
//		CPPUNIT_ASSERT_EQUAL( expected_y, y);
		CPPUNIT_ASSERT( (expected_y - y).norm() < BASSOTOLERANCE);
	}
}

void BackprojectionMatrixUnitTest::adjointTest()
{
	const unsigned int num_angles = 4;
	const unsigned int num_offsets = 5;
	const unsigned int num_pixels = 5;
	Eigen::VectorXd f(num_pixels*num_pixels);
	f.setRandom();
	Eigen::VectorXd g(num_angles*num_offsets);
	g.setRandom();
//	DiscretizedBackprojection RadonMatrix(num_pixels, num_pixels, num_angles, num_offsets);
	BackprojectionMatrix backprojection(num_pixels, num_pixels, num_angles, num_offsets);
//	// testing for adjoint correctness
	{
		const Eigen::MatrixXd& Rs = backprojection.getMatrix();
		const Eigen::MatrixXd R = backprojection.getMatrix().transpose();
		const Eigen::VectorXd Rf = R*f;
		const Eigen::VectorXd Rsg = Rs*g;

//		std::cout << "Got f=" << f.transpose() << std::endl;
//		std::cout << "Got Rf=" << Rf.transpose() << std::endl;
//		std::cout << "Got g=" << g.transpose() << std::endl;
//		std::cout << "Got Rsg=" << Rsg.transpose() << std::endl;
//		std::cout << "Then, <Rf,g> = " << g.dot(Rf)
//				<< " and <f,Rsg> = " << f.dot(Rsg) << std::endl;

		CPPUNIT_ASSERT( fabs(g.dot(Rf) - f.dot(Rsg)) < BASSOTOLERANCE );
	}

//	{
//
//		const Eigen::MatrixXd& R = RadonMatrix.getMatrix();
//		const Eigen::MatrixXd Rs = backprojection.getMatrix();
//		const Eigen::VectorXd Rf = R*f;
//		const Eigen::VectorXd Rsg = Rs*g;
//
//		std::cout << "Got f=" << f.transpose() << std::endl;
//		std::cout << "Got Rf=" << Rf.transpose() << std::endl;
//		std::cout << "Got g=" << g.transpose() << std::endl;
//		std::cout << "Got Rsg=" << Rsg.transpose() << std::endl;
//		std::cout << "Then, <Rf,g> = " << g.dot(Rf)
//				<< " and <f,Rsg> = " << f.dot(Rsg) << std::endl;
//
//		std::cout << "R=" << R << std::endl;
//		std::cout << "Rs'=" << Rs.transpose() << std::endl;
//		std::cout << R-Rs.transpose() << std::endl;
//
//		CPPUNIT_ASSERT( fabs(g.dot(Rf) - f.dot(Rsg)) < BASSOTOLERANCE );
//	}
}
