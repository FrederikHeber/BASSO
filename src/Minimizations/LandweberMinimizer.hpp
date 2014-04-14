/*
 * LandweberMinimizer.hpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#ifndef LANDWEBERMINIMIZER_HPP_
#define LANDWEBERMINIMIZER_HPP_

#include <Eigen/Dense>

#include "Minimizations/GeneralMinimizer.hpp"

#include "DualityMapping.hpp"
#include "LpNorm.hpp"
#include "SmoothnessModulus.hpp"

class LandweberMinimizer : public GeneralMinimizer
{
public:
	LandweberMinimizer(
			const double _NormX,
			const double _NormY,
			const double _PowerY,
			const double _Delta,
			const double _C,
			const unsigned int _maxiter
			);
	~LandweberMinimizer() {}

	GeneralMinimizer::ReturnValues operator()(
			const Eigen::VectorXd &_x0,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y
			) const;

private:
	/** Calculates tau for modulus of smoothness such that modulus
	 * over tau matches given \a _lambda
	 */
	double calculateMatchingTau(
			const SmoothnessModulus &_modul,
			const double _lambda) const;

public:
	/** Calculate residual \a _A * \a _x0 - \a _y in given norm \a _NormY.
	 *
	 * \param _x0 current iteration point
	 * \param _A matrix of inverse problem
	 * \param _y right-hand side
	 * \param _residual residual vector, updated after call
	 * \return norm of residual
	 */
	double calculateResidual(
			const Eigen::VectorXd &_x0,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y,
			Eigen::VectorXd &_residual
			) const;

	/** Calculate optimal step width via a line search.
	 *
	 * i.e. along u we look for alpha to minimize residual
	 *
	 * @param _x current position
	 * @param _u direction for line search
	 * @param _A matrix
	 * @param _y right-hand side
	 * @param _alpha initial guess for step width
	 * @return optimal step \f$ \alpha \f$ width to minimize
	 * 		\f$ || A (x - \alpha u) - y ||_Y \f$
	 */
	double calculateOptimalStepwidth(
			const Eigen::VectorXd &_x,
			const Eigen::VectorXd &_u,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y,
			const double _alpha = 0) const;

public:
	//!> Lp norm of space X: p
	const double val_NormX;
	//!> Lp norm of space Y: r
	const double val_NormY;
	//!> Lp norm of dual space to X: q
	const double val_DualNormX;
	//!> power of dual map J_p
	const double PowerX;
	//!> power of dual map J_r
	const double PowerY;
	//!> power of dual map J_q
	const double DualPowerX;
	//!> magnitude of noise
	const double Delta;
	//!> maximum number of iterations in outer loop
	const unsigned int MaxOuterIterations;
	//!> tolerance for objects in space X
	const double TolX;
	//!> tolerance for objects in space Y
	const double TolY;
	//!> tolerance for Fun
	const double TolFun;
	//!> positive dampening constant for iteration
	const double C;
	//!> norm object for space X
	const LpNorm NormX;
	//!> norm object for space Y
	const LpNorm NormY;
	//!> norm object for dual space X^*
	const LpNorm DualNormX;
	//!> duality mapping object for space X
	const DualityMapping J_p;
	//!> duality mapping object for dual space X^*
	const DualityMapping J_q;
	//!> duality mapping object for space Y (single-valued)
	const DualityMapping j_r;
	//!> smoothness modulus object for dual Space X^*
	const SmoothnessModulus modul;
};


#endif /* LANDWEBERMINIMIZER_HPP_ */
