/*
 * SequentialSubspaceMinimizerNoise.hpp
 *
 *  Created on: Jun 03, 2014
 *      Author: heber
 */

#ifndef SEQUENTIALSUBSPACEMINIMIZERNOISE_HPP_
#define SEQUENTIALSUBSPACEMINIMIZERNOISE_HPP_

#include <Eigen/Dense>

#include "Minimizations/GeneralMinimizer.hpp"

class Database;

/** This class implements the sequential subspace optimization by [Schöpfer,
 * Schuster,Louis, 2006].
 *
 */
class SequentialSubspaceMinimizerNoise : public GeneralMinimizer
{
public:
	SequentialSubspaceMinimizerNoise(
			const double _NormX,
			const double _NormY,
			const double _PowerX,
			const double _PowerY,
			const double _Delta,
			const unsigned int _maxiter,
			const Eigen::VectorXd &_solution,
			Database &_database,
			const unsigned int _outputsteps
			);
	~SequentialSubspaceMinimizerNoise() {}

	/** Setter for tau.
	 *
	 * This is to have a definite place where tau is changed. Hence,
	 * it is const and cannot accidentally be changed in the code, but
	 * it can still be set after the instance has been created.
	 *
	 * @param _tau new value of tau, tau in (1, infty)
	 */
	void setTau(const double _tau);

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


#endif /* SEQUENTIALSUBSPACEMINIMIZERNOISE_HPP_ */
