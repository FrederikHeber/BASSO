/*
 * LandweberFixedStepWidth.hpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */

#ifndef LANDWEBERFIXEDSTEPWIDTH_HPP_
#define LANDWEBERSTEPFIXEDWIDTH_HPP_

#include "BassoConfig.h"

#include "Minimizations/Functions/SmoothnessFunctional.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Minimizers/StepWidths/DetermineStepWidth.hpp"
#include "Minimizations/types.hpp"

/** Implements fixed landweber step width.
 *
 */
struct LandweberFixedStepWidth : public DetermineStepWidth
{
	LandweberFixedStepWidth(
			const InverseProblem_ptr_t &_problem,
			const double _C);

	const double operator()(
			const SpaceElement_ptr_t &_dualx,
			const SpaceElement_ptr_t &_u,
			const SpaceElement_ptr_t &_solution,
			const SpaceElement_ptr_t &_residual,
			const double _residuum,
			const double _TolX,
			const double _alpha) const;


private:
	double calculateMatchingTau(const double _lambda) const;

private:
	const InverseProblem_ptr_t problem;
	//!> smoothness modulus object for dual Space X^*
	const SmoothnessModulus modul;

	const double C;
	double G;
	double modulus_at_one;
	double ANorm;
};


#endif /* LANDWEBERFIXEDSTEPWIDTH_HPP_ */
