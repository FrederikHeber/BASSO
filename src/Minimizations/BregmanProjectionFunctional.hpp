/*
 * BregmanProjectionFunctional.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#ifndef BREGMANPROJECTIONFUNCTIONAL_HPP_
#define BREGMANPROJECTIONFUNCTIONAL_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "MatrixIO/OperationCounter.hpp"
#include "Minimizations/DualityMapping.hpp"
#include "Minimizations/LpNorm.hpp"

/** Functor to calculate BregmanProjectionFunctional functional/distance.
 *
 */
class BregmanProjectionFunctional
{
public:
	/** Constructor for BregmanProjectionFunctional.
	 *
	 * @param _lpdualnorm norm object of dual space
	 * @param _J_q duality mapping from dual space to space
	 */
	BregmanProjectionFunctional(
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
				> &_ScalarVectorProduct);
	~BregmanProjectionFunctional() {}

	/** Implements BregmanProjectionFunctional functional.
	 *
	 * \param _t coefficients to column vectors of search space
	 * \param _x current dual of solution
	 * \param _U search space spanned by column vectors
	 * \param _alpha offsets of affine subspace
	 * \param _q power of the weight of the duality mapping
	 */
	double operator()(
			const Eigen::VectorXd &_t,
			const Eigen::VectorXd &_dualx,
			const Eigen::MatrixXd &_U,
			const Eigen::VectorXd &_alpha,
			const double _q
			);

	/** Implements BregmanProjectionFunctional functional.
	 *
	 * \param _t coefficients to column vectors of search space
	 * \param _x current dual of solution
	 * \param _U search space spanned by column vectors
	 * \param _alpha offsets of affine subspace
	 * \param _q power of the weight of the duality mapping
	 */
	Eigen::VectorXd gradient(
			const Eigen::VectorXd &_t,
			const Eigen::VectorXd &_dualx,
			const Eigen::MatrixXd &_U,
			const Eigen::VectorXd &_alpha,
			const double _q
			);

private:
	//!> lp Norm object
	const LpNorm &lpdualnorm;
	//!> DualityMapping object
	const DualityMapping &J_q;
	//!> counting and timing object for MatrixVectorMultiplication
	const OperationCounter<
			const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type,
			const Eigen::MatrixBase<Eigen::MatrixXd>&,
			const Eigen::MatrixBase<Eigen::VectorXd>&
			> &MatrixVectorProduct;
	//!> counting and timing object for VectorVectorMultiplication
	const OperationCounter<
			Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
			const Eigen::MatrixBase<Eigen::VectorXd>&,
			const Eigen::MatrixBase<Eigen::VectorXd>&
			> &ScalarVectorProduct;
};


#endif /* BREGMANPROJECTIONFUNCTIONAL_HPP_ */
