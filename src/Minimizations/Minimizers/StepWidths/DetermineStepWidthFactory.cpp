/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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
		const ResidualFunctional::calculateResidual_t &_residualizer,
		const Mapping_ptr_t &_J_q)
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
			new DynamicRegularizedL1NormStepWidth(_problem, _J_q)
			);
		break;
	}
	return instance;
}


