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
			const StoppingCriterion::ptr_t &_stopping_criteria,
			const unsigned int _outputsteps,
			const LastNSearchDirections::OrthogonalizationType _orthogonalization_type
			);

	~SequentialSubspaceMinimizerNoise() {}

	GeneralMinimizer::ReturnValues operator()(
			const InverseProblem_ptr_t &_problem,
			const SpaceElement_ptr_t &_startvalue,
			const SpaceElement_ptr_t &_dualstartvalue,
			const SpaceElement_ptr_t &_truesolution
			);
};


#endif /* SEQUENTIALSUBSPACEMINIMIZERNOISE_HPP_ */
