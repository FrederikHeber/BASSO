/*
 * SequentialSubspaceMinimizer.hpp
 *
 *  Created on: Mar 28, 2014
 *      Author: heber
 */

#ifndef SEQUENTIALSUBSPACEMINIMIZER_HPP_
#define SEQUENTIALSUBSPACEMINIMIZER_HPP_

#include <Eigen/Dense>

/** This class implements the sequential subspace optimization by [Schöpfer,
 * Schuster,Louis, 2006].
 *
 */
class SequentialSubspaceMinimizer
{
public:
	SequentialSubspaceMinimizer();
	~SequentialSubspaceMinimizer() {}

	/** Internal structure for return values.
	 *
	 */
	struct ReturnValues
	{
		//!> solution vector
		Eigen::VectorXd solution;
		//!> remaining residuum
		double residuum;
		//!> number of outer iterations till solution
		int NumberOuterIterations;
	};

	ReturnValues operator()(
			const Eigen::VectorXd &_x0,
			const unsigned int _NormX,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y,
			const unsigned int _NormY,
			const double _PowerY,
			const double _Delta
			);

private:
	// internal variables

//	//!> number of columns (M)
//	unsigned int NoCols;
//	//!> number of rows (N)
//	unsigned int NoRows;

	// constants

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
