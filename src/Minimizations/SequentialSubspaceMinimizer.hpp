/*
 * SequentialSubspaceMinimizer.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef SEQUENTIALSUBSPACEMINIMIZER_HPP_
#define SEQUENTIALSUBSPACEMINIMIZER_HPP_

#include <Eigen/Dense>

#include "Minimizations/GeneralMinimizer.hpp"

/** This class implements the sequential subspace optimization by [Schöpfer,
 * Schuster,Louis, 2006].
 *
 */
class SequentialSubspaceMinimizer : public GeneralMinimizer
{
public:
	SequentialSubspaceMinimizer(
			const unsigned int _NormX,
			const unsigned int _NormY,
			const double _PowerY,
			const double _Delta
			);
	~SequentialSubspaceMinimizer() {}

	GeneralMinimizer::ReturnValues operator()(
			const Eigen::VectorXd &_x0,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y
			) const;

private:
	// internal variables

//	//!> number of columns (M)
//	unsigned int NoCols;
//	//!> number of rows (N)
//	unsigned int NoRows;

	// constants

	//!> Lp norm of space X: p
	const double val_NormX;
	//!> Lp norm of space Y: r
	const double val_NormY;
	//!> power of dual map J_r
	const double PowerY;
	//!> magnitude of noise
	const double Delta;
	//!> maximum number of iterations in outer loop
	const unsigned int MaxOuterIterations;
	//!> tolerance for x
	const double TolX;
	//!> tolerance for Fun
	const double TolFun;
	//!> regularization parameter for discrepancy principle, tau > 1
	const double tau;
};


#endif /* SEQUENTIALSUBSPACEMINIMIZER_HPP_ */
