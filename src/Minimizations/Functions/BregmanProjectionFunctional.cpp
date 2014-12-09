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
#include <numeric>

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
		const std::vector<double> &_t,
		const SpaceElement_ptr_t &_dualx,
		const std::vector<SpaceElement_ptr_t> &_U,
		const std::vector<double> &_alpha
		) const
{
	assert ( _t.size() == _U.size() );
	assert ( _t.size() == _alpha.size() );
	// x=x-U*t;
	SpaceElement_ptr_t resx = _dualx->getSpace()->createElement();
	*resx = std::inner_product(_t.begin(), _t.end(),_U.begin(),resx);
	*resx = _dualx - resx;
	// fval=1/q*norm(x,p)^q+alpha'*t;
	double alpha_times_t = 0.;
	alpha_times_t = std::inner_product(
			_t.begin(), _t.end(),
			_alpha.begin(),
			alpha_times_t);
	const double fval =
			1./dualpower * ::pow(dualnorm(resx), dualpower)
			+ alpha_times_t;
	return fval;
}

std::vector<double> BregmanProjectionFunctional::gradient(
		const std::vector<double> &_t,
		const SpaceElement_ptr_t &_dualx,
		const std::vector<SpaceElement_ptr_t> &_U,
		const std::vector<double> &_alpha
		) const
{
	assert ( _t.size() == _U.size() );
	assert ( _t.size() == _alpha.size() );
	// x=x-U*t;
	SpaceElement_ptr_t resx = _dualx->getSpace()->createElement();
	*resx = std::inner_product(_t.begin(), _t.end(),_U.begin(),resx);
	*resx = _dualx - resx;
	std::vector<double> gval(_alpha);
	const SpaceElement_ptr_t dual_resx = J_q(resx);
	for (size_t i=0;i<_t.size();++i)
		gval[i] -= _U[i] * dual_resx;
	return gval;
}
