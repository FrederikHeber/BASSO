/*
 * SequentialSubspaceMinimizerNoise.hpp
 *
 *  Created on: Jun 03, 2014
 *      Author: heber
 */

#ifndef SEQUENTIALSUBSPACEMINIMIZERNOISE_HPP_
#define SEQUENTIALSUBSPACEMINIMIZERNOISE_HPP_

#include "BassoConfig.h"

#include "Minimizations/Minimizers/SequentialSubspaceMinimizer.hpp"

class Database;

/** This class implements the sequential subspace optimization by [Schöpfer,
 * Schuster,Louis, 2006].
 *
 */
class SequentialSubspaceMinimizerNoise : public SequentialSubspaceMinimizer
{
public:
	SequentialSubspaceMinimizerNoise(
			const InverseProblem_ptr_t &_inverseproblem,
			const double _Delta,
			const unsigned int _maxiter,
			const unsigned int _maxinneriter,
			Database &_database,
			const unsigned int _outputsteps,
			const bool _orthogonal_directions
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
			const InverseProblem_ptr_t &_problem,
			const SpaceElement_ptr_t &_startvalue,
			const SpaceElement_ptr_t &_dualstartvalue,
			const SpaceElement_ptr_t &_truesolution
			);

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
