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

#include "Minimizations/BregmanProjectionFunctional.hpp"
#include "Minimizations/Functions/MinimizationFunctional.hpp"

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
	BregmanProjectionFunctional &bregman;
	const Eigen::VectorXd &x;
	const Eigen::MatrixXd &U;
	const Eigen::VectorXd &alpha;
	const double q;

	/** Constructor to initialize refs.
	 *
	 */
	HyperplaneProjection(
		BregmanProjectionFunctional &_bregman,
		const Eigen::VectorXd &_x,
		const Eigen::MatrixXd &_U,
		const Eigen::VectorXd &_alpha,
		const double _q
		);

	double operator()(const Eigen::VectorXd &_value) const;

	const Eigen::VectorXd gradient(const Eigen::VectorXd &_value) const;

	void convertToInternalType(
			Eigen::VectorXd &_t,
			const gsl_vector * const x) const;

	void convertFromInternalType(
			const Eigen::VectorXd &_t,
			gsl_vector * x) const;
};



#endif /* HYPERPLANEPROJECTION_HPP_ */
