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
	BregmanFunctional(const double _p, const double _tolerance = BASSOTOLERANCE);
	~BregmanFunctional() {}

	/** Implements BregmanFunctional functional.
	 *
	 * \param _t
	 * \param _x vector
	 * \param _U
	 * \param _alpha
	 * \param _q power of the weight of the duality mapping
	 */
	double operator()(
			const Eigen::VectorXd &_t,
			const Eigen::VectorXd &_x,
			const Eigen::MatrixXd &_U,
			const Eigen::VectorXd &_alpha,
			const unsigned int _q
			);

	/** Implements BregmanFunctional functional.
	 *
	 * \param _t
	 * \param _x vector
	 * \param _U
	 * \param _alpha
	 * \param _q power of the weight of the duality mapping
	 */
	Eigen::VectorXd gradient(
			const Eigen::VectorXd &_t,
			const Eigen::VectorXd &_x,
			const Eigen::MatrixXd &_U,
			const Eigen::VectorXd &_alpha,
			const unsigned int _q
			);

private:
	//!> value p of the Lp norm
	const double p;
	//!> lp Norm object
	LpNorm lpnorm;
	//!> DualityMapping object
	DualityMapping J_p;
};


#endif /* BREGMANFUNCTIONAL_HPP_ */
