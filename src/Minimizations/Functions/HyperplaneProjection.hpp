/*
 * HyperplaneProjection.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef HYPERPLANEPROJECTION_HPP_
#define HYPERPLANEPROJECTION_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Functions/BregmanProjectionFunctional.hpp"
#include "Minimizations/Functions/MinimizationFunctional.hpp"
#include "Minimizations/types.hpp"

/** Structure containing all parameters to call BregmanProjectionFunctional functions.
 *
 * This is required to use function minimization that only allows
 * to pass a void* pointer to pass on information to the function to be
 * minimized.
 *
 * It calculates the projection onto the intersection of hyperplanes.
 *
 * \sa BregmanProjectionFunctional
 *
 */
struct HyperplaneProjection :
		public MinimizationFunctional<Eigen::VectorXd>
{
	typedef typename MinimizationFunctional<Eigen::VectorXd>::array_type array_type;

	BregmanProjectionFunctional &bregman;
	const Eigen::VectorXd &x;
	const Eigen::MatrixXd &U;
	const Eigen::VectorXd &alpha;

	/** Constructor to initialize refs.
	 *
	 */
	HyperplaneProjection(
		BregmanProjectionFunctional &_bregman,
		const Eigen::VectorXd &_x,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha
		);

	/** Constructor to initialize refs.
	 *
	 */
	HyperplaneProjection(
		BregmanProjectionFunctional &_bregman,
		const SpaceElement_ptr_t &_x,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha
		);

	double function(const Eigen::VectorXd &_value) const;

	const Eigen::VectorXd gradient(const Eigen::VectorXd &_value) const;

	void convertInternalTypeToArrayType(
			const Eigen::VectorXd &_t,
			array_type & _x
			) const;

	void convertArrayTypeToInternalType(
			const array_type & _x,
			Eigen::VectorXd &_t
			) const;
};



#endif /* HYPERPLANEPROJECTION_HPP_ */
