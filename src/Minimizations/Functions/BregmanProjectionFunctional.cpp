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

#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/LpDualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"

BregmanProjectionFunctional::BregmanProjectionFunctional(
		const Norm &_dualnorm,
		const PowerTypeDualityMapping &_J_q,
		const double _dualpower,
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
	dualpower(_dualpower),
	dualnorm(_dualnorm),
	J_q(_J_q),
	MatrixVectorProduct(_MatrixVectorProduct),
	ScalarVectorProduct(_ScalarVectorProduct)
{}

BregmanProjectionFunctional::BregmanProjectionFunctional(
		const InverseProblem_ptr_t &_problem,
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
	dualpower(_problem->x->getSpace()->getDualSpace()->getDualityMapping()->getPower()),
	dualnorm(*_problem->x->getSpace()->getDualSpace()->getNorm()),
	J_q(dynamic_cast<const PowerTypeDualityMapping&>(
			*_problem->x->getSpace()->getDualSpace()->getDualityMapping())
			),
	MatrixVectorProduct(_MatrixVectorProduct),
	ScalarVectorProduct(_ScalarVectorProduct)
{}

double BregmanProjectionFunctional::operator()(
		const Eigen::VectorXd &_t,
		const Eigen::VectorXd &_dualx,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha
		) const
{
	// x=x-U*t;
	const Eigen::VectorXd resx = _dualx - MatrixVectorProduct(_U,_t);
	// fval=1/q*norm(x,p)^q+alpha'*t;
	const Eigen::VectorXd alpha_transposed = _alpha.transpose();
	const double fval =
			1./dualpower * ::pow(dualnorm(resx), dualpower)
			+ ScalarVectorProduct(alpha_transposed, _t);
	return fval;
}

double BregmanProjectionFunctional::operator()(
		const Eigen::VectorXd &_t,
		const SpaceElement_ptr_t &_dualx,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha
		) const
{
	// x=x-U*t;
	const Eigen::VectorXd resx =
			_dualx->getVectorRepresentation() - MatrixVectorProduct(_U,_t);
	// fval=1/q*norm(x,p)^q+alpha'*t;
	const Eigen::VectorXd alpha_transposed = _alpha.transpose();
	const double fval =
			1./dualpower * ::pow(dualnorm(resx), dualpower)
			+ ScalarVectorProduct(alpha_transposed, _t);
	return fval;
}

Eigen::VectorXd BregmanProjectionFunctional::gradient(
		const Eigen::VectorXd &_t,
		const Eigen::VectorXd &_dualx,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha
		) const
{
	const Eigen::VectorXd resx = _dualx - MatrixVectorProduct(_U, _t);
	const Eigen::MatrixXd &U_transposed = _U.transpose();
	const Eigen::VectorXd gval =
			_alpha -
			MatrixVectorProduct(U_transposed, J_q(resx));

	return gval;
}

Eigen::VectorXd BregmanProjectionFunctional::gradient(
		const Eigen::VectorXd &_t,
		const SpaceElement_ptr_t &_dualx,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha
		) const
{
	const Eigen::VectorXd resx =
			_dualx->getVectorRepresentation() - MatrixVectorProduct(_U, _t);
	const Eigen::MatrixXd &U_transposed = _U.transpose();
	const Eigen::VectorXd gval =
			_alpha -
			MatrixVectorProduct(U_transposed, J_q(resx));

	return gval;
}
