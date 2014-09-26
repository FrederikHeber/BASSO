/*
 * VectorProjectionUnitTest.cpp
 *
 *  Created on: Sep 10, 2014
 *      Author: heber
 */

#include "VectorProjectionUnitTest.hpp"

#include <boost/bind.hpp>
#include <Eigen/Dense>
#include <iostream>

#include "Log/Logging.hpp"
#include "Minimizations/Functions/VectorProjection.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( VectorProjectionUnitTest );

// static instances
const double VectorProjectionUnitTest::tolerance = 1e-8;

void VectorProjectionUnitTest::setUp()
{
	// VectorProjection uses logging
	logging::core::get()->set_filter
			(
					logging::trivial::severity >= logging::trivial::info
			);
	startLogging();

	scalar_vector_fctor =
			boost::bind(
					static_cast<Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType
						(Eigen::MatrixBase<Eigen::VectorXd>::*)(const Eigen::MatrixBase<Eigen::VectorXd>&) const>(
								&Eigen::MatrixBase<Eigen::VectorXd>::dot),
								_1, _2
			);
	ScalarVectorProduct = new OperationCounter<
			Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
			const Eigen::MatrixBase<Eigen::VectorXd>&,
			const Eigen::MatrixBase<Eigen::VectorXd>&
			>(scalar_vector_fctor);
}


void VectorProjectionUnitTest::tearDown()
{
	delete ScalarVectorProduct;
}

void VectorProjectionUnitTest::oneoneNorm()
{
	const unsigned int dim = 5;
	const double p=1.1;
	const NormedSpace_ptr_t SpaceX =
			NormedSpaceFactory::createLpInstance(dim, p, p);
	const PowerTypeDualityMapping &Jp =
			dynamic_cast<const PowerTypeDualityMapping &>(
					*SpaceX->getDualityMapping());
	const Norm &NormX = *SpaceX->getNorm();
	VectorProjection projector(NormX, Jp, p, *ScalarVectorProduct);
	Eigen::VectorXd projectonto(dim);
	Eigen::VectorXd tobeprojected(dim);

	// parallel vectors, unit length
	{
		projectonto << 1.5,0,0,0,0;
		tobeprojected << 1.,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, more than unit length
	for (double length = .1; length < 2.; length+=0.01) {
//		std::cout << "Current length is " << length << std::endl;
		projectonto << 1.5,0,0,0,0;
		tobeprojected << length,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e3);
	}

	// parallel vectors, less than unit length
	{
		projectonto << .1,0,0,0,0;
		tobeprojected << .2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, mixed
	{
		projectonto << 1.1,0,0,0,0;
		tobeprojected << .2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, other way mixed
	{
		projectonto << .27,0,0,0,0;
		tobeprojected << 2.2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, unit length
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected << 1,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// orthogonal vectors, unit length
	{
		projectonto << 1.5,0,0,0,0;
		tobeprojected << 0,1,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( tmp.second/NormX(tobeprojected) < tolerance*1e2);
	}

	// orthogonal vectors
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected << 0,0,1,1,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( tmp.second/NormX(tobeprojected) < tolerance*1e2);
	}

	// slanted vectors
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
		}
	}

	// slanted vectors, smaller power
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
		}
	}

	// slanted vectors, larger power
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
		}
	}
}

void VectorProjectionUnitTest::onefiveNorm()
{
	const unsigned int dim = 5;
	const double p=1.5;
	const NormedSpace_ptr_t SpaceX =
			NormedSpaceFactory::createLpInstance(dim, p, p);
	const PowerTypeDualityMapping &Jp =
			dynamic_cast<const PowerTypeDualityMapping &>(
					*SpaceX->getDualityMapping());
	const Norm &NormX = *SpaceX->getNorm();
	VectorProjection projector(NormX, Jp, p, *ScalarVectorProduct);
	Eigen::VectorXd projectonto(dim);
	Eigen::VectorXd tobeprojected(dim);

	// parallel vectors, unit length
	{
		projectonto << 1.5,0,0,0,0;
		tobeprojected << 1.,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, more than unit length
	for (double length = .1; length < 2.; length+=0.01) {
//		std::cout << "Current length is " << length << std::endl;
		projectonto << 1.5,0,0,0,0;
		tobeprojected << length,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, less than unit length
	{
		projectonto << .1,0,0,0,0;
		tobeprojected << .2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, mixed
	{
		projectonto << 1.1,0,0,0,0;
		tobeprojected << .2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, other way mixed
	{
		projectonto << .27,0,0,0,0;
		tobeprojected << 2.2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, unit length
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected << 1,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// orthogonal vectors, unit length
	{
		projectonto << 1.5,0,0,0,0;
		tobeprojected << 0,1,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( tmp.second/NormX(tobeprojected) < tolerance*1e2);
	}

	// orthogonal vectors
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected << 0,0,1,1,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( tmp.second/NormX(tobeprojected) < tolerance*1e2);
	}

	// slanted vectors
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
		}
	}

	// slanted vectors, smaller power
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
		}
	}

	// slanted vectors, larger power
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
		}
	}
}

void VectorProjectionUnitTest::twoNorm()
{
	const unsigned int dim = 5;
	const double p=2.;
	const NormedSpace_ptr_t SpaceX =
			NormedSpaceFactory::createLpInstance(dim, p, p);
	const PowerTypeDualityMapping &Jp =
			dynamic_cast<const PowerTypeDualityMapping &>(
					*SpaceX->getDualityMapping());
	const Norm &NormX = *SpaceX->getNorm();
	VectorProjection projector(NormX, Jp, p, *ScalarVectorProduct);
	Eigen::VectorXd projectonto(dim);
	Eigen::VectorXd tobeprojected(dim);

	// parallel vectors, unit length
	{
		projectonto << 1.5,0,0,0,0;
		tobeprojected << 1.,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
		const double angle = scalar_vector_fctor(projectonto, tobeprojected)
				/ NormX(projectonto) / NormX(tobeprojected);
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - angle) < tolerance);
	}

	// parallel vectors, more than unit length
	for (double length = .1; length < 2.; length+=0.01) {
//		std::cout << "Current length is " << length << std::endl;
		projectonto << 1.5,0,0,0,0;
		tobeprojected << length,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
		const double angle = scalar_vector_fctor(projectonto, tobeprojected)
				/ NormX(projectonto) / NormX(tobeprojected);
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - angle) < tolerance);
	}

	// parallel vectors, less than unit length
	{
		projectonto << .1,0,0,0,0;
		tobeprojected << .2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
		const double angle = scalar_vector_fctor(projectonto, tobeprojected)
				/ NormX(projectonto) / NormX(tobeprojected);
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - angle) < tolerance);
	}

	// parallel vectors, mixed
	{
		projectonto << 1.1,0,0,0,0;
		tobeprojected << .2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
		const double angle = scalar_vector_fctor(projectonto, tobeprojected)
				/ NormX(projectonto) / NormX(tobeprojected);
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - angle) < tolerance);
	}

	// parallel vectors, other way mixed
	{
		projectonto << .27,0,0,0,0;
		tobeprojected << 2.2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
		const double angle = scalar_vector_fctor(projectonto, tobeprojected)
				/ NormX(projectonto) / NormX(tobeprojected);
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - angle) < tolerance);
	}

	// parallel vectors, unit length
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected << 1,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
		const double angle = scalar_vector_fctor(projectonto, tobeprojected)
				/ NormX(projectonto) / NormX(tobeprojected);
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - angle) < tolerance);
	}

	// orthogonal vectors, unit length
	{
		projectonto << 1.5,0,0,0,0;
		tobeprojected << 0,1,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( tmp.second/NormX(tobeprojected) < tolerance*1e2);
		const double angle = scalar_vector_fctor(projectonto, tobeprojected)
				/ NormX(projectonto) / NormX(tobeprojected);
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - angle) < tolerance);
	}

	// orthogonal vectors
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected << 0,0,1,1,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( tmp.second/NormX(tobeprojected) < tolerance*1e2);
		const double angle = scalar_vector_fctor(projectonto, tobeprojected)
				/ NormX(projectonto) / NormX(tobeprojected);
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - angle) < tolerance);
	}

	// slanted vectors
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		std::vector<double> euclidianangles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
			euclidianangles[i] = scalar_vector_fctor(projectonto, tobeprojected)
					/ NormX(projectonto) / NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		CPPUNIT_ASSERT( fabs( euclidianangles[0] - angles[0]) < tolerance);
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
			CPPUNIT_ASSERT( fabs( euclidianangles[i] - angles[i]) < tolerance);
		}
	}

	// slanted vectors, smaller power
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		std::vector<double> euclidianangles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
			euclidianangles[i] = scalar_vector_fctor(projectonto, tobeprojected)
					/ NormX(projectonto) / NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		CPPUNIT_ASSERT( fabs( euclidianangles[0] - angles[0]) < tolerance);
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
			CPPUNIT_ASSERT( fabs( euclidianangles[i] - angles[i]) < tolerance);
		}
	}

	// slanted vectors, larger power
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		std::vector<double> euclidianangles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
			euclidianangles[i] = scalar_vector_fctor(projectonto, tobeprojected)
					/ NormX(projectonto) / NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		CPPUNIT_ASSERT( fabs( euclidianangles[0] - angles[0]) < tolerance);
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
			CPPUNIT_ASSERT( fabs( euclidianangles[i] - angles[i]) < tolerance);
		}
	}
}

void VectorProjectionUnitTest::sixNorm()
{
	const unsigned int dim = 5;
	const double p=6.;
	const NormedSpace_ptr_t SpaceX =
			NormedSpaceFactory::createLpInstance(dim, p, p);
	const PowerTypeDualityMapping &Jp =
			dynamic_cast<const PowerTypeDualityMapping &>(
					*SpaceX->getDualityMapping());
	const Norm &NormX = *SpaceX->getNorm();
	VectorProjection projector(NormX, Jp, p, *ScalarVectorProduct);
	Eigen::VectorXd projectonto(dim);
	Eigen::VectorXd tobeprojected(dim);

	// parallel vectors, unit length
	{
		projectonto << 1.5,0,0,0,0;
		tobeprojected << 1.,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, more than unit length
	for (double length = .1; length < 2.; length+=0.01) {
//		std::cout << "Current length is " << length << std::endl;
		projectonto << 1.5,0,0,0,0;
		tobeprojected << length,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
		std::cout << "Minimum is " << tmp.first
				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, less than unit length
	{
		projectonto << .1,0,0,0,0;
		tobeprojected << .2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, mixed
	{
		projectonto << 1.1,0,0,0,0;
		tobeprojected << .2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, other way mixed
	{
		projectonto << .27,0,0,0,0;
		tobeprojected << 2.2,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// parallel vectors, unit length
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected << 1,0,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( fabs( tmp.second/NormX(tobeprojected) - 1.) < tolerance*1e2);
	}

	// orthogonal vectors, unit length
	{
		projectonto << 1.5,0,0,0,0;
		tobeprojected << 0,1,0,0,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( tmp.second/NormX(tobeprojected) < tolerance*1e2);
	}

	// orthogonal vectors
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected << 0,0,1,1,0;
		// minimum, minimizer
		const std::pair<double, double> tmp = projector(
				projectonto,
				tobeprojected,
				tolerance
				);
//		std::cout << "Minimum is " << tmp.first
//				<< ", with minimizer " << tmp.second << std::endl;
		CPPUNIT_ASSERT( tmp.second/NormX(tobeprojected) < tolerance*1e2);
	}

	// slanted vectors
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
		}
	}

	// slanted vectors, smaller power
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
		}
	}

	// slanted vectors, larger power
	{
		projectonto << 1.,0,0,0,0;
		tobeprojected = Eigen::VectorXd::Zero(dim);
		std::vector<double> angles(dim, 0.);
		for (unsigned int i=0; i<dim; ++i) {
			tobeprojected[i] = 1.;
			// minimum, minimizer
			const std::pair<double, double> tmp = projector(
					projectonto,
					tobeprojected,
					tolerance
					);
//			std::cout << "Minimum is " << tmp.first
//					<< ", with minimizer " << tmp.second << std::endl;
			angles[i] = tmp.second/NormX(tobeprojected);
		}
		// make sure that angle is monotonically decreasing
		for (unsigned int i=1; i<dim; ++i) {
//			std::cout << "Comparing " << angles[i-1] << " with "
//					<< angles[i] << std::endl;
			CPPUNIT_ASSERT( angles[i-1] > angles[i]);
		}
	}
}