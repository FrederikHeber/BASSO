/*
 * NemirovskyDirection.cpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NemirovskyDirection.hpp"

#include <cassert>

#include "Log/Logging.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"

NemirovskyDirection::NemirovskyDirection(
		const NormedSpace_ptr_t &_SearchDirectionSpace_ptr,
		const OperationCounter<
			const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type,
			const Eigen::MatrixBase<Eigen::MatrixXd>&,
			const Eigen::MatrixBase<Eigen::VectorXd>&
			> &_MatrixVectorProduct_subspace,
		const OperationCounter<
			Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
			const Eigen::MatrixBase<Eigen::VectorXd>&,
			const Eigen::MatrixBase<Eigen::VectorXd>&
			> &_ScalarVectorProduct_subspace) :
		Searchspace(_SearchDirectionSpace_ptr, 2, _ScalarVectorProduct_subspace),
		MatrixVectorProduct_subspace(_MatrixVectorProduct_subspace),
		ScalarVectorProduct_subspace(_ScalarVectorProduct_subspace)
{}

void NemirovskyDirection::update(
		const SpaceElement_ptr_t &_newdir,
		const double _alpha,
		const SpaceElement_ptr_t &_dual_iterate,
		const SpaceElement_ptr_t &_iterate
		)
{
	assert( _newdir->getSpace() == SearchDirectionSpace_ptr );
	assert( _dual_iterate->getSpace() == SearchDirectionSpace_ptr );

	// update search direction
	U.col(0) = _newdir->getVectorRepresentation();
	alphas(0) = _alpha;

	// update Nemirovsky direction
	U.col(1) = _dual_iterate->getVectorRepresentation();
	alphas(1) = ScalarVectorProduct_subspace(
			_dual_iterate->getVectorRepresentation(),
			_iterate->getVectorRepresentation());
}
