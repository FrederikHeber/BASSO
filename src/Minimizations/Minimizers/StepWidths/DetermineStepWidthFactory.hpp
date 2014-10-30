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
	static DetermineStepWidth_ptr_t createInstance(
			const InverseProblem_ptr_t &_problem,
			const bool _useOptimalStepWidth,
			const double _C,
			const ResidualFunctional::calculateResidual_t &_residualizer);
};


#endif /* DETERMINESTEPWIDTHFACTORY_HPP_ */
