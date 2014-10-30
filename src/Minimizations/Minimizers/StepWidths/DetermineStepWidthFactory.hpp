/*
 * DetermineStepWidthFactory.hpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */

#ifndef DETERMINESTEPWIDTHFACTORY_HPP_
#define DETERMINESTEPWIDTHFACTORY_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

#include "Minimizations/Functions/ResidualFunctional.hpp"

/** Factory (function) for step width calculation.
 *
 */
struct DetermineStepWidthFactory
{
	//!> enumeration of all available step width calculators
	enum stepwidth_enumeration {
		LandweberFixes=0,
		MinimizingResidual=1,
		ConstantRegularizedL1Norm=2,
		DynamicRegularizedL1Norm=3,
		MAX_stepwidth_enumeration
	};

	static DetermineStepWidth_ptr_t createInstance(
			const InverseProblem_ptr_t &_problem,
			const enum stepwidth_enumeration _stepwidth_type,
			const double _C,
			const ResidualFunctional::calculateResidual_t &_residualizer);
};


#endif /* DETERMINESTEPWIDTHFACTORY_HPP_ */
