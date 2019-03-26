/*
 * ResidualMinimizingStepwidth.hpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */

#ifndef RESIDUALMINIMIZINGSTEPWIDTH_HPP_
#define RESIDUALMINIMIZINGSTEPWIDTH_HPP_

#include "BassoConfig.h"

#include "Minimizations/Minimizers/StepWidths/DetermineStepWidth.hpp"
#include "Minimizations/Functions/ResidualFunctional.hpp"
#include "Minimizations/types.hpp"

/** This implements a default step size where new residuum is minimized
 * via a Brent algorithm (does not require derivative).
 */
struct ResidualMinimizingStepwidth : public DetermineStepWidth
{
	ResidualMinimizingStepwidth(
			const InverseProblem_ptr_t &_problem,
			const ResidualFunctional::calculateResidual_t &_residualizer
			) :
		problem(_problem),
		residualizer(_residualizer)
	{}

	virtual ~ResidualMinimizingStepwidth() {}

	const double operator()(
			const SpaceElement_ptr_t &_dualx,
			const SpaceElement_ptr_t &_u,
			const SpaceElement_ptr_t &_solution,
			const SpaceElement_ptr_t &_residual,
			const double _residuum,
			const double _TolX,
			const double _alpha) const;

private:

	const InverseProblem_ptr_t problem;
	const ResidualFunctional::calculateResidual_t &residualizer;
};



#endif /* RESIDUALMINIMIZINGSTEPWIDTH_HPP_ */
