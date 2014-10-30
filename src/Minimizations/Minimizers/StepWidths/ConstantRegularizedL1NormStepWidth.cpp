/*
 * ConstantRegularizedL1NormStepWidth.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "ConstantRegularizedL1NormStepWidth.hpp"

#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"

ConstantRegularizedL1NormStepWidth::ConstantRegularizedL1NormStepWidth(
		const InverseProblem_ptr_t &_problem) :
		ANorm_reciprocal_sqr(
				1./::pow(
						dynamic_cast<const LinearMapping &>(*_problem->A).Norm(),
						2))
{}


