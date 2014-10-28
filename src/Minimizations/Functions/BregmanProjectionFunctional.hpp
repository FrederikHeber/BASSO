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
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"

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
				> &_ScalarVectorProduct);
	~BregmanProjectionFunctional() {}

	/** Implements BregmanProjectionFunctional functional.
	 *
	 * \param _t coefficients to column vectors of search space
	 * \param _x current dual of solution
	 * \param _U search space spanned by column vectors
	 * \param _alpha offsets of affine subspace
	 */
	double operator()(
			const Eigen::VectorXd &_t,
			const Eigen::VectorXd &_dualx,
			const Eigen::MatrixXd &_U,
			const Eigen::VectorXd &_alpha
			) const;

	/** Implements BregmanProjectionFunctional functional.
	 *
	 * \param _t coefficients to column vectors of search space
	 * \param _x current dual of solution
	 * \param _U search space spanned by column vectors
	 * \param _alpha offsets of affine subspace
	 */
	Eigen::VectorXd gradient(
			const Eigen::VectorXd &_t,
			const Eigen::VectorXd &_dualx,
			const Eigen::MatrixXd &_U,
			const Eigen::VectorXd &_alpha
			) const;

private:
	//!> power type of the weight function of \a J_q
	const double dualpower;
	//!> lp Norm object
	const Norm &dualnorm;
	//!> LpDualityMapping object
	const PowerTypeDualityMapping &J_q;
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
