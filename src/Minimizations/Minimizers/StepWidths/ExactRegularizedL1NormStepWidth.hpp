/*
 * ExactRegularizedL1NormStepWidth.hpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */

#ifndef EXACTREGULARIZEDL1NORMSTEPWIDTH_HPP_
#define EXACTREGULARIZEDL1NORMSTEPWIDTH_HPP_

#include "BassoConfig.h"

#include "Minimizations/Minimizers/StepWidths/DetermineStepWidth.hpp"

struct ExactRegularizedL1NormStepWidth :
		public DetermineStepWidth
{
	ExactRegularizedL1NormStepWidth(
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
};



#endif /* EXACTREGULARIZEDL1NORMSTEPWIDTH_HPP_ */
