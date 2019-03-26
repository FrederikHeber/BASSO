/*
 * DiscretizedRadonMatrixUnitTest.cpp
 *
 *  Created on: Jul 07, 2015
 *      Author: heber
 */

#include "BassoConfig.h"

#include "ComputerTomography/DiscretizedRadon/unittests/DiscretizedRadonMatrixUnitTest.hpp"

#include <Eigen/Dense>

#include "Log/Logging.hpp"

#include <ComputerTomography/DiscretizedRadon/BackprojectionMatrix.hpp>
#include "ComputerTomography/DiscretizedRadon/DiscretizedRadonMatrix.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( DiscretizedRadonMatrixUnitTest );

void DiscretizedRadonMatrixUnitTest::setUp()
{
	// DiscretizedRadonMatrix uses logging
	boost::log::core::get()->set_filter
			(
					boost::log::trivial::severity >= boost::log::trivial::debug
			);
	startLogging();
}


void DiscretizedRadonMatrixUnitTest::tearDown()
{
}

const double domain_length=(1.-(-1.));

void DiscretizedRadonMatrixUnitTest::singlePixel_few_angles_few_offsetsTest()
{
	// constants
	const int num_pixel_x = 1;
	const int num_pixel_y = 1;
	const int num_pixels = num_pixel_x*num_pixel_y;

	const int num_offsets = 1;
	const int num_angles = 2;
	const int num_measurements = num_angles*num_offsets;

	// construct matrix
	DiscretizedRadonMatrix RadonMatrix(
			num_pixel_x,
			num_pixel_y,
			num_angles,
			num_offsets);
	const Eigen::MatrixXd &A = RadonMatrix.getMatrix();

	// single pixel with white value one
	{
		// create expected value
		Eigen::VectorXd x(num_pixels);
		x.setZero();
		x(0) = 1.;

		Eigen::VectorXd expected_y(num_measurements);
		expected_y.setZero();
		expected_y(0) = 1.*domain_length;
		expected_y(1) = 1.*domain_length;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}
	// single pixel with grey value one
	{
		// create expected value
		Eigen::VectorXd x(num_pixels);
		x.setZero();
		x(0) = .5;

		Eigen::VectorXd expected_y(num_measurements);
		expected_y.setZero();
		expected_y(0) = .5*domain_length;
		expected_y(1) = .5*domain_length;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}
}

inline double IntersectionLengthThroughPixelAtAngle(
		const double _angle)
{
	const double radian = _angle/180*M_PI;
	return sqrt(1.+::pow(tan(radian),2));
}

void DiscretizedRadonMatrixUnitTest::singlePixel_many_angles_few_offsetsTest()
{
	// constants
	const int num_pixel_x = 1;
	const int num_pixel_y = 1;
	const int num_pixels = num_pixel_x*num_pixel_y;

	const int num_offsets = 1;
	const int num_angles = 8;
	const int num_measurements = num_angles*num_offsets;

	// construct matrix
	DiscretizedRadonMatrix RadonMatrix(
			num_pixel_x,
			num_pixel_y,
			num_angles,
			num_offsets);
	const Eigen::MatrixXd &A = RadonMatrix.getMatrix();

	const double fourth_pi_quarter_distance =
			domain_length*IntersectionLengthThroughPixelAtAngle(22.5);
	const double half_pi_quarter_distance =
			domain_length*IntersectionLengthThroughPixelAtAngle(45.);

	// single pixel with white value one
	{
		// create expected value
		Eigen::VectorXd x(num_pixels);
		x.setZero();
		x(0) = 1.;

		Eigen::VectorXd expected_y(num_measurements);
		expected_y.setZero();
		expected_y(0) = 1.*domain_length;
		expected_y(1) = fourth_pi_quarter_distance;
		expected_y(2) = half_pi_quarter_distance;
		expected_y(3) = fourth_pi_quarter_distance;
		expected_y(4) = 1.*domain_length;
		expected_y(5) = fourth_pi_quarter_distance;
		expected_y(6) = half_pi_quarter_distance;
		expected_y(7) = fourth_pi_quarter_distance;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}
	// single pixel with grey value one
	{
		// create expected value
		Eigen::VectorXd x(num_pixels);
		x.setZero();
		x(0) = .5;

		Eigen::VectorXd expected_y(num_measurements);
		expected_y.setZero();
		expected_y(0) = .5*domain_length;
		expected_y(1) = .5*fourth_pi_quarter_distance;
		expected_y(2) = .5*half_pi_quarter_distance;
		expected_y(3) = .5*fourth_pi_quarter_distance;
		expected_y(4) = .5*domain_length;
		expected_y(5) = .5*fourth_pi_quarter_distance;
		expected_y(6) = .5*half_pi_quarter_distance;
		expected_y(7) = .5*fourth_pi_quarter_distance;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}
}

void DiscretizedRadonMatrixUnitTest::singlePixel_few_angles_many_offsetsTest()
{
	// constants
	const int num_pixel_x = 1;
	const int num_pixel_y = 1;
	const int num_pixels = num_pixel_x*num_pixel_y;

	const int num_offsets = 3;
	const int num_angles = 2;
	const int num_measurements = num_angles*num_offsets;

	// construct matrix
	DiscretizedRadonMatrix RadonMatrix(
			num_pixel_x,
			num_pixel_y,
			num_angles,
			num_offsets);
	const Eigen::MatrixXd &A = RadonMatrix.getMatrix();

	// single pixel with white value one
	{
		// create expected value
		Eigen::VectorXd x(num_pixels);
		x.setZero();
		x(0) = 1.;

		Eigen::VectorXd expected_y(num_measurements);
		expected_y.setZero();
		expected_y(0) = 1.*domain_length;
		expected_y(1) = 1.*domain_length;
		expected_y(2) = 1.*domain_length;
		expected_y(3) = 1.*domain_length;
		expected_y(4) = 1.*domain_length;
		expected_y(5) = 1.*domain_length;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}
	// single pixel with grey value one
	{
		// create expected value
		Eigen::VectorXd x(num_pixels);
		x.setZero();
		x(0) = .5;

		Eigen::VectorXd expected_y(num_measurements);
		expected_y.setZero();
		expected_y(0) = .5*domain_length;
		expected_y(1) = .5*domain_length;
		expected_y(2) = .5*domain_length;
		expected_y(3) = .5*domain_length;
		expected_y(4) = .5*domain_length;
		expected_y(5) = .5*domain_length;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}
}

void DiscretizedRadonMatrixUnitTest::ninePixel_few_angles_few_offsetsTest()
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
	DiscretizedRadonMatrix RadonMatrix(
			num_pixel_x,
			num_pixel_y,
			num_angles,
			num_offsets);
	const Eigen::MatrixXd &A = RadonMatrix.getMatrix();

	// nine pixel with single white pixel in center
	{
		// create expected value
		Eigen::VectorXd x(num_pixels);
		x.setZero();
		x(1 + (1 * num_pixel_y)) = 1.;

		// single ray, intensity is basically a constant for every angle
		Eigen::VectorXd expected_y(num_measurements);
		expected_y.setZero();
		expected_y(0) = 1./3.*domain_length;
		expected_y(1) = 1./3.*domain_length;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}
}


void DiscretizedRadonMatrixUnitTest::ninePixel_few_angles_many_offsetsTest()
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
	DiscretizedRadonMatrix RadonMatrix(
			num_pixel_x,
			num_pixel_y,
			num_angles,
			num_offsets);
	const Eigen::MatrixXd &A = RadonMatrix.getMatrix();

	// nine pixel with single white pixel in center
	{
		// create expected value
		Eigen::VectorXd x(num_pixels);
		x.setZero();
		x(1 + (1 * num_pixel_y)) = 1.;

		// single ray does not hit central pixel at all
		Eigen::VectorXd expected_y(num_measurements);
		expected_y.setZero();
		expected_y(0) = 0.*domain_length;
		expected_y(1) = 1./3.*domain_length;
		expected_y(2) = 0.*domain_length;
		expected_y(3) = 0.*domain_length;
		expected_y(4) = 1./3.*domain_length;
		expected_y(5) = 0.*domain_length;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}

	// nine pixel with single white central horizontal line
	{
		// create expected value
		Eigen::VectorXd x(num_pixels);
		x.setZero();
		for (unsigned int i=0;i<num_pixel_x;++i)
			x(i + (1 * num_pixel_y)) = 1.;

		Eigen::VectorXd expected_y(num_measurements);
		expected_y.setZero();
		expected_y(0) = 0.*domain_length;
		expected_y(1) = 1.*domain_length;
		expected_y(2) = 0.*domain_length;
		expected_y(3) = 1./3.*domain_length;
		expected_y(4) = 1./3.*domain_length;
		expected_y(5) = 1./3.*domain_length;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}

	// nine pixel with single white central vertical line
	{
		// create expected value
		Eigen::VectorXd x(num_pixels);
		x.setZero();
		for (unsigned int i=0;i<num_pixel_x;++i)
			x(1 + (i * num_pixel_y)) = 1.;

		Eigen::VectorXd expected_y(num_measurements);
		expected_y.setZero();
		expected_y(0) = 1./3.*domain_length;
		expected_y(1) = 1./3.*domain_length;
		expected_y(2) = 1./3.*domain_length;
		expected_y(3) = 0.*domain_length;
		expected_y(4) = 1.*domain_length;
		expected_y(5) = 0.*domain_length;

		// calculate y
		Eigen::VectorXd y = A * x;
		CPPUNIT_ASSERT_EQUAL( expected_y, y );
	}
}

void DiscretizedRadonMatrixUnitTest::adjointTest()
{
	const unsigned int num_angles = 4;
	const unsigned int num_offsets = 5;
	const unsigned int num_pixels = 5;
	Eigen::VectorXd f(num_pixels*num_pixels);
	f.setRandom();
	Eigen::VectorXd g(num_angles*num_offsets);
	g.setRandom();
	DiscretizedRadonMatrix RadonMatrix(num_pixels, num_pixels, num_angles, num_offsets);
//	BackprojectionMatrix backprojection(num_pixels, num_pixels, num_angles, num_offsets);
	// testing for adjoint correctness
	{
		const Eigen::MatrixXd& R = RadonMatrix.getMatrix();
		const Eigen::MatrixXd Rs = RadonMatrix.getMatrix().transpose();
		const Eigen::VectorXd Rf = R*f;
		const Eigen::VectorXd Rsg = Rs*g;

//		std::cout << "Got f=" << f.transpose() << std::endl;
//		std::cout << "Got Rf=" << Rf.transpose() << std::endl;
//		std::cout << "Got g=" << g.transpose() << std::endl;
//		std::cout << "Got Rsg=" << Rsg.transpose() << std::endl;
//		std::cout << "Then, <Rf,g> = " << g.dot(Rf)
//				<< " and <f,Rsg> = " << f.dot(Rsg) << std::endl;

		// The following should work only for true adjoint
//		CPPUNIT_ASSERT( fabs(g.dot(Rf) - f.dot(Rsg)) < BASSOTOLERANCE );
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
