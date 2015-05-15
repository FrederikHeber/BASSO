/*
 * DynamicRegularizedL1NormStepWidth.hpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */

#ifndef DYNAMICREGULARIZEDL1NORMSTEPWIDTH_HPP_
#define DYNAMICREGULARIZEDL1NORMSTEPWIDTH_HPP_

#include "BassoConfig.h"

#include "Minimizations/Minimizers/StepWidths/DetermineStepWidth.hpp"
#include "Minimizations/types.hpp"

/** This implements the "dynamic" stepsize according to
 * [Lorenz et al., '13, §3.1]
 */
struct DynamicRegularizedL1NormStepWidth :
		public DetermineStepWidth
{
	DynamicRegularizedL1NormStepWidth(
			const InverseProblem_ptr_t &_problem);

	const double operator()(
			const SpaceElement_ptr_t &_dualx,
			const SpaceElement_ptr_t &_u,
			const SpaceElement_ptr_t &_solution,
			const SpaceElement_ptr_t &_residual,
			const double _residuum,
			const double _TolX,
			const double _alpha
			) const;

private:
	const InverseProblem_ptr_t problem;
	const Mapping_ptr_t A_adjoint;
};


#endif /* DYNAMICREGULARIZEDL1NORMSTEPWIDTH_HPP_ */
