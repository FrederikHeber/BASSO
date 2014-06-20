/*
 * BregmanFunctional.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#ifndef BREGMANFUNCTIONAL_HPP_
#define BREGMANFUNCTIONAL_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/DualityMapping.hpp"
#include "Minimizations/LpNorm.hpp"

/** Functor to calculate BregmanFunctional functional/distance.
 *
 */
class BregmanFunctional
{
public:
	/** Constructor for BregmanFunctional.
	 *
	 * @param _lpdualnorm norm object of dual space
	 * @param _J_q duality mapping from dual space to space
	 */
	BregmanFunctional(
			const LpNorm &_lpdualnorm,
			const DualityMapping &_J_q);
	~BregmanFunctional() {}

	/** Implements BregmanFunctional functional.
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

	/** Implements BregmanFunctional functional.
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
};


#endif /* BREGMANFUNCTIONAL_HPP_ */
