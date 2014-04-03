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
