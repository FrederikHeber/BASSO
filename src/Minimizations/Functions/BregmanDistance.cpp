/*
 * BregmanDistance.cpp
 *
 *  Created on: Apr 25, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "BregmanDistance.hpp"

#include <cmath>
#include <Eigen/Dense>
#include <limits>

#include "Log/Logging.hpp"
#include "Math/Helpers.hpp"
#include "Minimizations/Mappings/LpDualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"

BregmanDistance::BregmanDistance(
		const Norm &_norm,
		const PowerTypeDualityMapping &_J_p,
		const double _power,
		const OperationCounter<
				Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
				const Eigen::MatrixBase<Eigen::VectorXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&
				>& _ScalarVectorProduct) :
			power(_power),
			norm(_norm),
			J_p(_J_p),
			ScalarVectorProduct(_ScalarVectorProduct)
{
	if ((power != 0.) && (power <= 1.))
		throw MinimizationIllegalValue_exception()
			<< MinimizationIllegalValue_name("power");
}

double BregmanDistance::operator()(
		const Eigen::VectorXd &_x,
		const Eigen::VectorXd &_y
		) const
{
	const Eigen::VectorXd dual_x = J_p(_x).transpose();
	return operator()(_x,_y, dual_x);
}

double BregmanDistance::operator()(
		const Eigen::VectorXd &_x,
		const Eigen::VectorXd &_y,
		const Eigen::VectorXd &_xdual
		) const
{

	BOOST_LOG_TRIVIAL(trace)
			<< "Calculating Bregman distance between "
			<< _x.transpose() << " and " << _y.transpose();
	double result = 0.;
	result += (1./Helpers::ConjugateValue(power)) * ::pow(norm(_x), power);
	result += (1./power) * ::pow(norm(_y), power);
	result -= ScalarVectorProduct(_xdual, _y);
	return result;
}
