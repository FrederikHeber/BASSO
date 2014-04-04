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

class LpNorm;
class SmoothnessModulus;

class LandweberMinimizer : public GeneralMinimizer
{
public:
	LandweberMinimizer();
	~LandweberMinimizer() {}

	GeneralMinimizer::ReturnValues operator()(
			const Eigen::VectorXd &_x0,
			const unsigned int _NormX,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y,
			const unsigned int _NormY,
			const double _PowerY,
			const double _Delta
			) const;

private:
	/** Calculates tau for modulus of smoothness such that modulus
	 * over tau matches given \a _lambda
	 */
	double calculateMatchingTau(
			SmoothnessModulus &_modul,
			const double _lambda) const;

	/** Calculate residual \a _A * \a _x0 - \a _y in given norm \a _NormY.
	 *
	 * \param _x0 current iteration point
	 * \param _A matrix of inverse problem
	 * \param _y right-hand side
	 * \param _NormY functor to calculate norm
	 * \param _residual residual vector, updated after call
	 * \return norm of residual
	 */
	double calculateResidual(
			const Eigen::VectorXd &_x0,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y,
			const LpNorm &_NormY,
			Eigen::VectorXd &_residual
			) const;

private:
	//!> maximum number of iterations in outer loop
	const unsigned int MaxOuterIterations;
	//!> tolerance for x
	const double TolX;
	//!> tolerance for Fun
	const double TolFun;
	//!> dampening constant for iteration
	const double C;
};


#endif /* LANDWEBERMINIMIZER_HPP_ */
