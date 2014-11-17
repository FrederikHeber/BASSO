/*
 * NemirovskyDirection.hpp
 *
 *  Created on: Nov 17, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_NEMIROVSKYDIRECTION_HPP_
#define MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_NEMIROVSKYDIRECTION_HPP_

#include "BassoConfig.h"

#include "Searchspace.hpp"

class NemirovskyDirection : public Searchspace
{
public:
	/** Constructor for class NemirovskyDirection.
	 *
	 * We have here two search directions, the current one
	 * and the current dual iterate.
	 *
	 * @param _SearchDirectionSpace_ptr search direction space (for checks)
	 * @param _MatrixVectorProduct_subspace counts for matrix-vector
	 * @param _ScalarVectorProduct_subspace counts for vector-vector
	 * products in the subspace (i.e. search space)
	 */
	NemirovskyDirection(
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
				> &_ScalarVectorProduct_subspace);

	/** This function performs the actual update of the search space.
	 *
	 * @param _newdir current search direction
	 * @param _alpha current alpha to this \a _newdir
	 * @param _dual_iterate current dual iterate
	 * @param _iterate current iterate
	 */
	void update(
			const SpaceElement_ptr_t &_newdir,
			const double _alpha,
			const SpaceElement_ptr_t &_dual_iterate,
			const SpaceElement_ptr_t &_iterate
			);

	/** Returns 0, as is always the search direction.
	 *
	 * @return 0
	 */
	const unsigned int getIndex() const
	{ return 0; }

private:
	//!> counter for the small matrix vector products in subspace
	const OperationCounter<
		const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type,
		const Eigen::MatrixBase<Eigen::MatrixXd>&,
		const Eigen::MatrixBase<Eigen::VectorXd>&
		> &MatrixVectorProduct_subspace;

	//!> counter for the small scalar products in subspace
	const OperationCounter<
		Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
		const Eigen::MatrixBase<Eigen::VectorXd>&,
		const Eigen::MatrixBase<Eigen::VectorXd>&
		> &ScalarVectorProduct_subspace;
};



#endif /* MINIMIZATIONS_MINIMIZERS_SEARCHSPACE_NEMIROVSKYDIRECTION_HPP_ */
