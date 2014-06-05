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
			const double _NormX,
			const double _NormY,
			const double _PowerX,
			const double _PowerY,
			const double _Delta,
			const unsigned int _maxiter,
			const Eigen::VectorXd &_solution,
			const unsigned int _outputsteps
			);
	~SequentialSubspaceMinimizer() {}

	/** Setter for tau.
	 *
	 * This is to have a definite place where tau is changed. Hence,
	 * it is const and cannot accidentally be changed in the code, but
	 * it can still be set after the instance has been created.
	 *
	 * @param _tau new value of tau, tau in (1, infty)
	 */
	void setTau(const double _tau);

	/** Setter for N, the number of search directions.
	 *
	 * This is to have a definite place where N is changed. Hence,
	 * it is const and cannot accidentally be changed in the code, but
	 * it can still be set after the instance has been created.
	 *
	 * @param _N new value of N, N in [1, infty)
	 */
	void setN(const unsigned int _N);

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

	//!> number of search directions
	const unsigned int N;
};


#endif /* SEQUENTIALSUBSPACEMINIMIZER_HPP_ */
