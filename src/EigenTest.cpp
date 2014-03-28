/*
 * EigenTest.cpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

#include "Minimizations/DualityMapping.hpp"

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
		DualityMapping<1> J_1(2);
		std::cout << "DualityMapping J_1 with weight 2 of v is ("
				<< J_1(v).transpose() << ")" << std::endl;
	}
	{
		DualityMapping<2> J_2(2);
		std::cout << "DualityMapping J_2 with weight 2 of v is ("
				<< J_2(v).transpose() << ")" << std::endl;
	}
	{
		DualityMapping<Eigen::Infinity> J_infty(2);
		std::cout << "DualityMapping J_infty with weight 2 of v is ("
				<< J_infty(v).transpose() << ")" << std::endl;
	}

}
