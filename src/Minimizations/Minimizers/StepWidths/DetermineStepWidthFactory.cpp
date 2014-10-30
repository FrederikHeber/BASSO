/*
 * DetermineStepWidthFactory.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "DetermineStepWidthFactory.hpp"

#include "Minimizations/Minimizers/StepWidths/ConstantRegularizedL1NormStepWidth.hpp"
#include "Minimizations/Minimizers/StepWidths/DynamicRegularizedL1NormStepWidth.hpp"
#include "Minimizations/Minimizers/StepWidths/LandweberFixedStepWidth.hpp"
#include "Minimizations/Minimizers/StepWidths/ResidualMinimizingStepwidth.hpp"

DetermineStepWidth_ptr_t
DetermineStepWidthFactory::createInstance(
		const InverseProblem_ptr_t &_problem,
		const enum stepwidth_enumeration _stepwidth_type,
		const double _C,
		const ResidualFunctional::calculateResidual_t &_residualizer)
{
	DetermineStepWidth_ptr_t instance;
	switch (_stepwidth_type) {
	default:
	case MinimizingResidual:
		instance = DetermineStepWidth_ptr_t(
				new ResidualMinimizingStepwidth(_problem,_residualizer)
				);
		break;
	case LandweberFixes:
		instance = DetermineStepWidth_ptr_t(
				new LandweberFixedStepWidth(_problem, _C)
				);
		break;
	case ConstantRegularizedL1Norm:
	instance = DetermineStepWidth_ptr_t(
			new ConstantRegularizedL1NormStepWidth(_problem)
			);
		break;
	case DynamicRegularizedL1Norm:
	instance = DetermineStepWidth_ptr_t(
			new DynamicRegularizedL1NormStepWidth(_problem)
			);
		break;
	}
	return instance;
}


