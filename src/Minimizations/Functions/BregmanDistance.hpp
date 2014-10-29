/*
 * BregmanDistance.hpp
 *
 *  Created on: Apr 25, 2014
 *      Author: heber
 */

#ifndef BREGMANDISTANCE_HPP_
#define BREGMANDISTANCE_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "MatrixIO/OperationCounter.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"

/** This implements a functor calculating the Bregman distance between
 * two points in a Lp space.
 *
 */
class BregmanDistance
{
public:
	BregmanDistance(
			const Norm &_norm,
			const PowerTypeDualityMapping &_J_p,
			const double _power,
			const OperationCounter<
					Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
					const Eigen::MatrixBase<Eigen::VectorXd>&,
					const Eigen::MatrixBase<Eigen::VectorXd>&
					>& _ScalarVectorProduct);
	~BregmanDistance() {}

	/** Calculate the Bregman distance between \a _x and \a _y.
	 *
	 * The Bregman distance reads as
	 * \f$ \Delta_p(x,y) = \tfrac 1 q ||x||^p + \tfrac 1 p ||y||^p - \langle J_{p}(x), y \rangle \f$
	 *
	 * Note that the argument \a p (and \a q with it) is ambigious because
	 * the lp-norm also has an argument p. For clarification, the p and q
	 * in the above formulas refer to the power type \a p of the duality
	 * mapping \f$ L_p \f$. Hence, q is conjugate to \a _power. We refer
	 * to  \f$ L_p \f$'s always as \a _power to distinguish clearly from
	 * the norm.
	 *
	 * That's also why \a p is set in the cstor where we create the norm
	 * object while the power type may be given in operator().
	 *
	 * @param _x first argument
	 * @param _y second argument
	 * @return Bregman distance between first and second argument
	 */
	double operator()(
			const SpaceElement_ptr_t &_x,
			const SpaceElement_ptr_t &_y
			) const
	{ return operator()(
			_x->getVectorRepresentation(),
			_y->getVectorRepresentation());
	}

	/** Calculate the Bregman distance between \a _x and \a _y.
	 *
	 * The Bregman distance reads as
	 * \f$ \Delta_p(x,y) = \tfrac 1 q ||x||^p + \tfrac 1 p ||y||^p - \langle J_{p}(x), y \rangle \f$
	 *
	 * Note that the argument \a p (and \a q with it) is ambigious because
	 * the lp-norm also has an argument p. For clarification, the p and q
	 * in the above formulas refer to the power type \a p of the duality
	 * mapping \f$ L_p \f$. Hence, q is conjugate to \a _power. We refer
	 * to  \f$ L_p \f$'s always as \a _power to distinguish clearly from
	 * the norm.
	 *
	 * That's also why \a p is set in the cstor where we create the norm
	 * object while the power type may be given in operator().
	 *
	 * @param _x first argument
	 * @param _y second argument
	 * @return Bregman distance between first and second argument
	 */
	double operator()(
			const Eigen::VectorXd &_x,
			const Eigen::VectorXd &_y
			) const;

	/** Calculate the Bregman distance between \a _x and \a _y with dual \a _xdual.
	 *
	 * The Bregman distance reads as
	 * \f$ \Delta_p(x,y) = \tfrac 1 q ||x||^p + \tfrac 1 p ||y||^p - \langle x^\ast, y \rangle \f$
	 *
	 * Note that the argument \a p (and \a q with it) is ambigious because
	 * the lp-norm also has an argument p. For clarification, the p and q
	 * in the above formulas refer to the power type \a p of the duality
	 * mapping \f$ L_p \f$. Hence, q is conjugate to \a _power. We refer
	 * to  \f$ L_p \f$'s always as \a _power to distinguish clearly from
	 * the norm.
	 *
	 * That's also why \a p is set in the cstor where we create the norm
	 * object while the power type may be given in operator().
	 *
	 * @param _x first argument
	 * @param _y second argument
	 * @param _xdual dual element to the first argument
	 * @return Bregman distance between first and second argument
	 */
	double operator()(
			const Eigen::VectorXd &_x,
			const Eigen::VectorXd &_y,
			const Eigen::VectorXd &_xdual
			) const;

	/** Calculate the Bregman distance between \a _x and \a _y with dual \a _xdual.
	 *
	 * The Bregman distance reads as
	 * \f$ \Delta_p(x,y) = \tfrac 1 q ||x||^p + \tfrac 1 p ||y||^p - \langle x^\ast, y \rangle \f$
	 *
	 * Note that the argument \a p (and \a q with it) is ambigious because
	 * the lp-norm also has an argument p. For clarification, the p and q
	 * in the above formulas refer to the power type \a p of the duality
	 * mapping \f$ L_p \f$. Hence, q is conjugate to \a _power. We refer
	 * to  \f$ L_p \f$'s always as \a _power to distinguish clearly from
	 * the norm.
	 *
	 * That's also why \a p is set in the cstor where we create the norm
	 * object while the power type may be given in operator().
	 *
	 * @param _x first argument
	 * @param _y second argument
	 * @param _xdual dual element to the first argument
	 * @return Bregman distance between first and second argument
	 */
	double operator()(
			const SpaceElement_ptr_t &_x,
			const SpaceElement_ptr_t &_y,
			const SpaceElement_ptr_t &_xdual
			) const
	{ return operator()(
			_x->getVectorRepresentation(),
			_y->getVectorRepresentation(),
			_xdual->getVectorRepresentation());
	}

private:
	//!> power type of the weight function of the duality mapping \a J_p
	const double power;
	//!> lp Norm object
	const Norm &norm;
	//!> LpDualityMapping object
	const PowerTypeDualityMapping &J_p;
	//!> counting and timing object for VectorVectorMultiplication
	const OperationCounter<
						Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
						const Eigen::MatrixBase<Eigen::VectorXd>&,
						const Eigen::MatrixBase<Eigen::VectorXd>&
						>& ScalarVectorProduct;
};


#endif /* BREGMANDISTANCE_HPP_ */
