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
#include <vector>

#include "MatrixIO/OperationCounter.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/types.hpp"

/** Functor to calculate BregmanProjectionFunctional functional/distance.
 *
 */
class BregmanProjectionFunctional
{
public:
	/** Constructor for BregmanProjectionFunctional.
	 *
	 * @param _dualnorm norm object of dual space
	 * @param _J_q duality mapping from dual space to space
	 * @param _dualpower power of this duality mapping
	 * @param _MatrixVectorProduct counter for matrix vector products
	 * @param _ScalarVectorProduct counter for scalar vector products
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

	/** Constructor for BregmanProjectionFunctional.
	 *
	 * @param _problem inverse problem with refs to norm and duality mapping
	 * @param _MatrixVectorProduct counter for matrix vector products
	 * @param _ScalarVectorProduct counter for scalar vector products
	 */
	BregmanProjectionFunctional(
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
			const std::vector<double> &_t,
			const SpaceElement_ptr_t &_dualx,
			const std::vector<SpaceElement_ptr_t> &_U,
			const std::vector<double> &_alpha
			) const;

	/** Implements BregmanProjectionFunctional functional.
	 *
	 * \param _t coefficients to column vectors of search space
	 * \param _x current dual of solution
	 * \param _U search space spanned by column vectors
	 * \param _alpha offsets of affine subspace
	 */
	std::vector<double> gradient(
			const std::vector<double> &_t,
			const SpaceElement_ptr_t &_dualx,
			const std::vector<SpaceElement_ptr_t> &_U,
			const std::vector<double> &_alpha
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
