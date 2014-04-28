/*
 * GeneralMinimizer.hpp
 *
 *  Created on: Apr 3, 2014
 *      Author: heber
 */

#ifndef GENERALMINIMIZER_HPP_
#define GENERALMINIMIZER_HPP_

#include <Eigen/Dense>

#include "DualityMapping.hpp"
#include "LpNorm.hpp"
#include "SmoothnessModulus.hpp"

/** This class describes the interface to a general minimizer.
 *
 */
class GeneralMinimizer
{
public:
	GeneralMinimizer(
			const double _NormX,
			const double _NormY,
			const double _PowerX,
			const double _PowerY,
			const double _Delta,
			const double _C,
			const unsigned int _maxiter,
			const unsigned int _outputsteps=0
			);
	virtual ~GeneralMinimizer() {}

	/** Internal structure for return values.
	 *
	 */
	struct ReturnValues
	{
		//!> solution vector
		Eigen::VectorXd solution;
		//!> residual vector
		Eigen::VectorXd residual;
		//!> remaining residuum, i.e. norm of residual
		double residuum;
		//!> number of outer iterations till solution
		int NumberOuterIterations;
	};

	virtual ReturnValues operator()(
			const Eigen::VectorXd &_x0,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y
			) const = 0;

protected:
	/** Internal helper function for specific Minimizers to print debugging
	 *  solutions.
	 *
	 * @param _solution intermediate solution
	 * @param _A forward matrix, required for printing projected solution
	 * @param _NumberOuterIterations current iteration count
	 */
	void printIntermediateSolution(
			const Eigen::VectorXd &_solution,
			const Eigen::MatrixXd &_A,
			unsigned int _NumberOuterIterations
			) const
;

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
	const int MaxOuterIterations;
	//!> tolerance for objects in space X
	const double TolX;
	//!> tolerance for objects in space Y
	const double TolY;
	//!> tolerance for Fun
	const double TolFun;
	//!> positive dampening constant for iteration
	const double C;
	//!> output solution each .. steps, 0 means never
	unsigned int outputsteps;

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
};


#endif /* GENERALMINIMIZER_HPP_ */
