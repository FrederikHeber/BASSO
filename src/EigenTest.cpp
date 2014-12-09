/*
 * EigenTest.cpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include <boost/function.hpp>

#include <boost/assign.hpp>
#include <boost/bind.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <numeric>
#include <sstream>

#include "MatrixIO/OperationCounter.hpp"
#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Mappings/Specifics/L1DualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LInfinityDualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LpDualityMapping.hpp"
#include "Minimizations/Norms/L1Norm.hpp"
#include "Minimizations/Norms/LInfinityNorm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Functions/BregmanProjectionFunctional.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

using namespace boost::assign;

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

	const double power = 2.;
	// testing of LpDualityMapping
	{
		const double p = 1.;
		NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						v.innerSize(), p, power);
		SpaceElement_ptr_t vElement =
				ElementCreator::create(SpaceX, v);
		Mapping_ptr_t J_1 = SpaceX->getDualityMapping();
		std::cout << "LpDualityMapping J_1 with weight 2 of v is ("
				<< (*J_1)(vElement) << ")" << std::endl;
	}
	{
		const double p = 2.;
		NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						v.innerSize(), p, power);
		SpaceElement_ptr_t vElement =
				ElementCreator::create(SpaceX, v);
		Mapping_ptr_t J_2 = SpaceX->getDualityMapping();
		std::cout << "LpDualityMapping J_2 with weight 2 of v is ("
				<< (*J_2)(vElement) << ")" << std::endl;
	}
	{
		const double p = std::numeric_limits<double>::infinity();
		NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						v.innerSize(), p, power);
		SpaceElement_ptr_t vElement =
				ElementCreator::create(SpaceX, v);
		Mapping_ptr_t J_infty = SpaceX->getDualityMapping();
		std::cout << "LpDualityMapping J_infty with weight 2 of v is ("
				<< (*J_infty)(vElement) << ")" << std::endl;
	}

	// testing of BregmanProjectionFunctional
	{
		const double p = 1.;
		const double power = 2.;
		NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						v.innerSize(), p, power);
		const Mapping_ptr_t &J_infty = SpaceX->getDualSpace()->getDualityMapping();
		BregmanProjectionFunctional bregman_1(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(*J_infty),
				J_infty->getPower());
		std::vector<double> t;
		t += 4,3;
		Eigen::VectorXd xtemp(2);
		xtemp << 4,3;
		SpaceElement_ptr_t x =
				ElementCreator::create(SpaceX, xtemp);
		Eigen::MatrixXd Utemp(2,2);
		Utemp << 1,0,0,1;
		std::vector<SpaceElement_ptr_t> U;
		for (size_t i=0;i<2;++i)
			U.push_back(
					ElementCreator::create(
							*SpaceX->getDualSpace(),
							Utemp.col(i)));
		assert( U.size() == (size_t)2 );
		std::vector<double> alpha;
		alpha += 1,0;
		const double fval = bregman_1(t,x,U,alpha);
		const std::vector<double> gval = bregman_1.gradient(t,x,U,alpha);
		std::stringstream gval_stream;
		std::copy(gval.begin(), gval.end(), std::ostream_iterator<double>(gval_stream, " "));
		std::cout << "BregmanProjectionFunctional bregman_2 of v is "
				<< fval << ","
				<< gval_stream.str() << "" << std::endl;
	}
	{
		const double p = 2.;
		const double power = 2.;
		NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						v.innerSize(), p, power);
		const Mapping_ptr_t &J_infty = SpaceX->getDualSpace()->getDualityMapping();
		BregmanProjectionFunctional bregman_2(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(*J_infty),
				J_infty->getPower());
		std::vector<double> t;
		t += 4,3;
		Eigen::VectorXd xtemp(2);
		xtemp << 4,3;
		SpaceElement_ptr_t x =
				ElementCreator::create(SpaceX, xtemp);
		Eigen::MatrixXd Utemp(2,2);
		Utemp << 1,0,0,1;
		std::vector<SpaceElement_ptr_t> U;
		for (size_t i=0;i<2;++i)
			U.push_back(
					ElementCreator::create(
							*SpaceX->getDualSpace(),
							Utemp.col(i)));
		assert( U.size() == (size_t)2 );
		std::vector<double> alpha;
		alpha += 1,0;
		const double fval = bregman_2(t,x,U,alpha);
		const std::vector<double> gval = bregman_2.gradient(t,x,U,alpha);
		std::stringstream gval_stream;
		std::copy(gval.begin(), gval.end(), std::ostream_iterator<double>(gval_stream, " "));
		std::cout << "BregmanProjectionFunctional bregman_2 of v is "
				<< fval << ","
				<< gval_stream.str() << "" << std::endl;
	}
	{
		const double p = std::numeric_limits<double>::infinity();
		const double power = 2.;
		NormedSpace_ptr_t SpaceX =
				NormedSpaceFactory::createLpInstance(
						v.innerSize(), p, power);
		const Mapping_ptr_t &J_infty = SpaceX->getDualSpace()->getDualityMapping();
		BregmanProjectionFunctional bregman_infty(
				*SpaceX->getDualSpace()->getNorm(),
				dynamic_cast<const PowerTypeDualityMapping &>(*J_infty),
				J_infty->getPower());
		std::vector<double> t;
		t += 4,3;
		Eigen::VectorXd xtemp(2);
		xtemp << 4,3;
		SpaceElement_ptr_t x =
				ElementCreator::create(SpaceX, xtemp);
		Eigen::MatrixXd Utemp(2,2);
		Utemp << 1,0,0,1;
		std::vector<SpaceElement_ptr_t> U;
		for (size_t i=0;i<2;++i)
			U.push_back(
					ElementCreator::create(
							*SpaceX->getDualSpace(),
							Utemp.col(i)));
		assert( U.size() == (size_t)2 );
		std::vector<double> alpha;
		alpha += 1,0;
		const double fval = bregman_infty(t,x,U,alpha);
		const std::vector<double> gval = bregman_infty.gradient(t,x,U,alpha);
		std::stringstream gval_stream;
		std::copy(gval.begin(), gval.end(), std::ostream_iterator<double>(gval_stream, " "));
		std::cout << "BregmanProjectionFunctional bregman_2 of v is "
				<< fval << ","
				<< gval_stream.str() << "" << std::endl;
	}
}
