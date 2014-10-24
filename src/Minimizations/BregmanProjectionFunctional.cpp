/*
 * BregmanProjectionFunctional.cpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#include "BassoConfig.h"

#include "BregmanProjectionFunctional.hpp"

#include <cmath>
#include <Eigen/Dense>

#include "Minimizations/DualityMappings/DualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"

BregmanProjectionFunctional::BregmanProjectionFunctional(
		const LpNorm &_lpdualnorm,
		const DualityMapping &_J_q,
		const OperationCounter<
			const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type,
			const Eigen::MatrixBase<Eigen::MatrixXd>&,
			const Eigen::MatrixBase<Eigen::VectorXd>&
			> &_MatrixVectorProduct,
		const OperationCounter<
			Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
			const Eigen::MatrixBase<Eigen::VectorXd>&,
			const Eigen::MatrixBase<Eigen::VectorXd>&
			> &_ScalarVectorProduct
		) :
	lpdualnorm(_lpdualnorm),
	J_q(_J_q),
	MatrixVectorProduct(_MatrixVectorProduct),
	ScalarVectorProduct(_ScalarVectorProduct)
{}

double BregmanProjectionFunctional::operator()(
		const Eigen::VectorXd &_t,
		const Eigen::VectorXd &_dualx,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha,
		const double _q
		)
{
	// x=x-U*t;
	const Eigen::VectorXd resx = _dualx - MatrixVectorProduct(_U,_t);
	// fval=1/q*norm(x,p)^q+alpha'*t;
	const Eigen::VectorXd alpha_transposed = _alpha.transpose();
	const double fval =
			1./(double)_q * ::pow(lpdualnorm(resx), _q)
			+ ScalarVectorProduct(alpha_transposed, _t);
	return fval;
}

Eigen::VectorXd BregmanProjectionFunctional::gradient(
		const Eigen::VectorXd &_t,
		const Eigen::VectorXd &_dualx,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha,
		const double _q
		)
{
	const Eigen::VectorXd resx = _dualx - MatrixVectorProduct(_U, _t);
	const Eigen::MatrixXd &U_transposed = _U.transpose();
	const Eigen::VectorXd gval =
			_alpha -
			MatrixVectorProduct(U_transposed, J_q(resx, _q));

	return gval;
}
