/*
 * DetermineStepWidthFactory.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "DetermineStepWidthFactory.hpp"

#include "Minimizations/Minimizers/StepWidths/LandweberFixedStepWidth.hpp"
#include "Minimizations/Minimizers/StepWidths/ResidualMinimizingStepwidth.hpp"

DetermineStepWidth_ptr_t
DetermineStepWidthFactory::createInstance(
		const InverseProblem_ptr_t &_problem,
		const bool _useOptimalStepWidth,
		const double _C,
		const ResidualFunctional::calculateResidual_t &_residualizer)
{
	DetermineStepWidth_ptr_t instance;
	if (_useOptimalStepWidth) {
		instance = DetermineStepWidth_ptr_t(
				new ResidualMinimizingStepwidth(_problem,_residualizer)
				);
	} else {
		instance = DetermineStepWidth_ptr_t(
				new LandweberFixedStepWidth(_problem, _C)
				);
	}
	return instance;
}


