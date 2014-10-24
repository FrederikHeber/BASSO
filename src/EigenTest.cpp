/*
 * EigenTest.cpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/function.hpp>

#include <boost/bind.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

#include "MatrixIO/OperationCounter.hpp"
#include "Minimizations/Mappings/L1DualityMapping.hpp"
#include "Minimizations/Mappings/LInfinityDualityMapping.hpp"
#include "Minimizations/Mappings/LpDualityMapping.hpp"
#include "Minimizations/Norms/L1Norm.hpp"
#include "Minimizations/Norms/LInfinityNorm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Functions/BregmanProjectionFunctional.hpp"

int main()
{
	Eigen::MatrixXd m(2,2);
	m(0,0) = 3;
	m(1,0) = 2.5;
	m(0,1) = -1;
	m(1,1) = m(1,0) + m(0,1);
	std::cout << "Here is the matrix m:\n" << m << std::endl;
	Eigen::VectorXd v(2);
	v(0) = 4;
	v(1) = v(0) - 1;
	std::cout << "Here is the vector v:\n" << v << std::endl;

	// testing Norm class
	std::cout << "Norms of v, L2: " << v.norm()
			<< ", L1: " << v.lpNorm<1>()
			<< ", l_infty: " << v.lpNorm<Eigen::Infinity>() << std::endl;

	// testing of LpDualityMapping
	{
		L1DualityMapping J_1;
		std::cout << "LpDualityMapping J_1 with weight 2 of v is ("
				<< J_1(v,2).transpose() << ")" << std::endl;
	}
	{
		LpDualityMapping J_2(2.);
		std::cout << "LpDualityMapping J_2 with weight 2 of v is ("
				<< J_2(v,2).transpose() << ")" << std::endl;
	}
	{
		LInfinityDualityMapping J_infty;
		std::cout << "LpDualityMapping J_infty with weight 2 of v is ("
				<< J_infty(v,2).transpose() << ")" << std::endl;
	}

	boost::function<
		const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type  (
				const Eigen::MatrixBase<Eigen::MatrixXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&
				)> matrix_vector_fctor =
			boost::bind(
					static_cast<const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type
						(Eigen::MatrixBase<Eigen::MatrixXd>::*)(const Eigen::MatrixBase<Eigen::VectorXd>&) const>(
								&Eigen::MatrixBase<Eigen::MatrixXd>::operator*),
								_1, _2
			)
	;
	const OperationCounter<
				const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type,
				const Eigen::MatrixBase<Eigen::MatrixXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&
				> MatrixVectorProduct(matrix_vector_fctor);
	boost::function<
			Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType (
					const Eigen::MatrixBase<Eigen::VectorXd>&,
					const Eigen::MatrixBase<Eigen::VectorXd>&)
					> scalar_vector_fctor =
			boost::bind(
					static_cast<Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType
						(Eigen::MatrixBase<Eigen::VectorXd>::*)(const Eigen::MatrixBase<Eigen::VectorXd>&) const>(
								&Eigen::MatrixBase<Eigen::VectorXd>::dot),
								_1, _2
			);
	const OperationCounter<
				Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
				const Eigen::MatrixBase<Eigen::VectorXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&
				> ScalarVectorProduct(scalar_vector_fctor);

	// testing of BregmanProjectionFunctional
	{
		L1Norm lpnorm;
		LInfinityNorm lpdualnorm;
		L1DualityMapping J_1;
		LInfinityDualityMapping J_infty;
		BregmanProjectionFunctional bregman_1(lpdualnorm, J_infty, MatrixVectorProduct, ScalarVectorProduct);
		Eigen::VectorXd t(2);
		t << 4,3;
		Eigen::VectorXd x(2);
		x << 4,3;
		Eigen::MatrixXd U(2,2);
		U << 1,0,0,1;
		Eigen::VectorXd alpha(2);
		alpha << 1,0;
		const unsigned int q = 2; 		// power of weight of duality mapping
		std::cout << "BregmanProjectionFunctional bregman_1 of v is "
				<< bregman_1(t,x,U,alpha,q) << ","
				<< bregman_1.gradient(t,x,U,alpha,q)<< "" << std::endl;
	}
	{
		LpNorm lpnorm(2);
		LpDualityMapping J_2(2);
		BregmanProjectionFunctional bregman_2(lpnorm, J_2, MatrixVectorProduct, ScalarVectorProduct);
		Eigen::VectorXd t(2);
		t << 4,3;
		Eigen::VectorXd x(2);
		x << 4,3;
		Eigen::MatrixXd U(2,2);
		U << 1,0,0,1;
		Eigen::VectorXd alpha(2);
		alpha << 1,0;
		const unsigned int q = 2; 		// power of weight of duality mapping
		std::cout << "BregmanProjectionFunctional bregman_2 of v is "
				<< bregman_2(t,x,U,alpha,q) << ","
				<< bregman_2.gradient(t,x,U,alpha,q) << "" << std::endl;
	}
	{
		LInfinityNorm lpnorm;
		L1Norm lpdualnorm;
		L1DualityMapping J_1;
		LInfinityDualityMapping J_infty;
		BregmanProjectionFunctional bregman_infty(lpdualnorm, J_1, MatrixVectorProduct, ScalarVectorProduct);
		Eigen::VectorXd t(2);
		t << 4,3;
		Eigen::VectorXd x(2);
		x << 4,3;
		Eigen::MatrixXd U(2,2);
		U << 1,0,0,1;
		Eigen::VectorXd alpha(2);
		alpha << 1,0;
		const unsigned int q = 2; 		// power of weight of duality mapping
		std::cout << "BregmanProjectionFunctional bregman_2 of v is "
				<< bregman_infty(t,x,U,alpha,q) << ","
				<< bregman_infty.gradient(t,x,U,alpha,q) << "" << std::endl;
	}
}
