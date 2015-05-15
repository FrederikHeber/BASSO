/*
 * ConstantRegularizedL1NormStepWidth.hpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */

#ifndef CONSTANTREGULARIZEDL1NORMSTEPWIDTH_HPP_
#define CONSTANTREGULARIZEDL1NORMSTEPWIDTH_HPP_

#include "BassoConfig.h"

#include "Minimizations/Minimizers/StepWidths/DetermineStepWidth.hpp"

/** This implements the "constant" stepsize according to
 * [Lorenz et al., '13, §3.1]
 */
struct ConstantRegularizedL1NormStepWidth :
		public DetermineStepWidth
{
	ConstantRegularizedL1NormStepWidth(
			const InverseProblem_ptr_t &_problem);

	const double operator()(
			const SpaceElement_ptr_t &_dualx,
			const SpaceElement_ptr_t &_u,
			const SpaceElement_ptr_t &_solution,
			const SpaceElement_ptr_t &_residual,
			const double _residuum,
			const double _TolX,
			const double _alpha
			) const
	{ return ANorm_reciprocal_sqr; }

private:
	const double ANorm_reciprocal_sqr;
};


#endif /* CONSTANTREGULARIZEDL1NORMSTEPWIDTH_HPP_ */
