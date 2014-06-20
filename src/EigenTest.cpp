/*
 * EigenTest.cpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

#include "Minimizations/DualityMapping.hpp"
#include "Minimizations/BregmanFunctional.hpp"

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

	// testing of DualityMapping
	{
		DualityMapping J_1(1);
		std::cout << "DualityMapping J_1 with weight 2 of v is ("
				<< J_1(v,2).transpose() << ")" << std::endl;
	}
	{
		DualityMapping J_2(2);
		std::cout << "DualityMapping J_2 with weight 2 of v is ("
				<< J_2(v,2).transpose() << ")" << std::endl;
	}
	{
		DualityMapping J_infty(LpNorm::Infinity);
		std::cout << "DualityMapping J_infty with weight 2 of v is ("
				<< J_infty(v,2).transpose() << ")" << std::endl;
	}

	// testing of BregmanFunctional
	{
		LpNorm lpnorm(1);
		LpNorm lpdualnorm(LpNorm::Infinity);
		DualityMapping J_1(1);
		DualityMapping J_infty(LpNorm::Infinity);
		BregmanFunctional bregman_1(lpdualnorm, J_infty);
		Eigen::VectorXd t(2);
		t << 4,3;
		Eigen::VectorXd x(2);
		x << 4,3;
		Eigen::MatrixXd U(2,2);
		U << 1,0,0,1;
		Eigen::VectorXd alpha(2);
		alpha << 1,0;
		const unsigned int q = 2; 		// power of weight of duality mapping
		std::cout << "BregmanFunctional bregman_1 of v is "
				<< bregman_1(t,x,U,alpha,q) << ","
				<< bregman_1.gradient(t,x,U,alpha,q)<< "" << std::endl;
	}
	{
		LpNorm lpnorm(2);
		DualityMapping J_2(2);
		BregmanFunctional bregman_2(lpnorm, J_2);
		Eigen::VectorXd t(2);
		t << 4,3;
		Eigen::VectorXd x(2);
		x << 4,3;
		Eigen::MatrixXd U(2,2);
		U << 1,0,0,1;
		Eigen::VectorXd alpha(2);
		alpha << 1,0;
		const unsigned int q = 2; 		// power of weight of duality mapping
		std::cout << "BregmanFunctional bregman_2 of v is "
				<< bregman_2(t,x,U,alpha,q) << ","
				<< bregman_2.gradient(t,x,U,alpha,q) << "" << std::endl;
	}
	{
		LpNorm lpnorm(LpNorm::Infinity);
		LpNorm lpdualnorm(1.);
		DualityMapping J_1(LpNorm::Infinity);
		DualityMapping J_infty(1.);
		BregmanFunctional bregman_infty(lpdualnorm, J_1);
		Eigen::VectorXd t(2);
		t << 4,3;
		Eigen::VectorXd x(2);
		x << 4,3;
		Eigen::MatrixXd U(2,2);
		U << 1,0,0,1;
		Eigen::VectorXd alpha(2);
		alpha << 1,0;
		const unsigned int q = 2; 		// power of weight of duality mapping
		std::cout << "BregmanFunctional bregman_2 of v is "
				<< bregman_infty(t,x,U,alpha,q) << ","
				<< bregman_infty.gradient(t,x,U,alpha,q) << "" << std::endl;
	}
}
