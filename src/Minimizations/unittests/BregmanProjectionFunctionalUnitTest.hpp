/*
 * BregmanProjectionFunctionalUnitTest.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef BregmanProjectionFunctionalUNITTEST_HPP_
#define BregmanProjectionFunctionalUNITTEST_HPP_

#include <cppunit/extensions/HelperMacros.h>

#include <Eigen/Dense>

#include <boost/function.hpp>

#include "MatrixIO/OperationCounter.hpp"

class BregmanProjectionFunctionalUnitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( BregmanProjectionFunctionalUnitTest );
  CPPUNIT_TEST( oneoneNorm );
  CPPUNIT_TEST( twoNorm );
  CPPUNIT_TEST( inftyNorm );
  CPPUNIT_TEST_SUITE_END();

public:
  void setUp();
  void tearDown();

  void oneoneNorm();
  void twoNorm();
  void inftyNorm();

private:
	boost::function<
		const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type  (
				const Eigen::MatrixBase<Eigen::MatrixXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&
				)> matrix_vector_fctor;
	OperationCounter<
		const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type,
		const Eigen::MatrixBase<Eigen::MatrixXd>&,
		const Eigen::MatrixBase<Eigen::VectorXd>&
		> *MatrixVectorProduct;

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

};


#endif /* BregmanProjectionFunctionalUNITTEST_HPP_ */
