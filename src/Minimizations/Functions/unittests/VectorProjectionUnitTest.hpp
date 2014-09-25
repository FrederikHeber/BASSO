/*
 * VectorProjectionUnitTest.hpp
 *
 *  Created on: Sep 10, 2014
 *      Author: heber
 */

#ifndef VECTORPROJECTIONUNITTEST_HPP_
#define VECTORPROJECTIONUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

#include <Eigen/Dense>

#include <boost/function.hpp>

#include "MatrixIO/OperationCounter.hpp"

class VectorProjectionUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( VectorProjectionUnitTest );
  CPPUNIT_TEST( oneoneNorm );
  CPPUNIT_TEST( onefiveNorm );
  CPPUNIT_TEST( twoNorm );
  CPPUNIT_TEST( sixNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneoneNorm();
  void onefiveNorm();
  void twoNorm();
  void sixNorm();

private:
	boost::function<
		Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType (
				const Eigen::MatrixBase<Eigen::VectorXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&)
				> scalar_vector_fctor;
	OperationCounter<
		Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
		const Eigen::MatrixBase<Eigen::VectorXd>&,
		const Eigen::MatrixBase<Eigen::VectorXd>&
		> *ScalarVectorProduct;

	//!> "global" tolerance threshold for minimization
	static const double tolerance;
};


#endif /* VECTORPROJECTIONUNITTEST_HPP_ */
