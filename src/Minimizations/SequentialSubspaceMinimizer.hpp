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
			const double _PowerX,
			const double _PowerY,
			const double _Delta,
			const double _C,
			const unsigned int _maxiter,
			const unsigned int _outputsteps
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

	//!> regularization parameter for discrepancy principle, tau > 1
	const double tau;
};


#endif /* SEQUENTIALSUBSPACEMINIMIZER_HPP_ */
