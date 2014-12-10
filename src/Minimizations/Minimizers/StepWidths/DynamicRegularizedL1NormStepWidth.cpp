/*
 * DynamicRegularizedL1NormStepWidth.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "DynamicRegularizedL1NormStepWidth.hpp"

#include <Eigen/Dense>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/NormFactory.hpp"

DynamicRegularizedL1NormStepWidth::DynamicRegularizedL1NormStepWidth(
		const InverseProblem_ptr_t &_problem) :
	problem(_problem),
	A_adjoint(_problem->A->getAdjointMapping())
{}

const double DynamicRegularizedL1NormStepWidth::operator()(
		const SpaceElement_ptr_t &_dualx,
		const SpaceElement_ptr_t &_u,
		const SpaceElement_ptr_t &_solution,
		const SpaceElement_ptr_t &_residual,
		const double _residuum,
		const double _TolX,
		const double _alpha
		) const
{
	// create L2 norm for measuring error
	const Norm_ptr_t l2norm = NormFactory::createLpInstance(
			problem->A->getTargetSpace(), 2.);

	const SpaceElement_ptr_t upper = problem->y->getSpace()->createElement();
	const LinearMapping &A = static_cast<const LinearMapping &>(*problem->A);
	*upper = A * _solution - problem->y;
	const SpaceElement_ptr_t lower = problem->y->getSpace()->createElement();
	*lower = A * upper;
	const double norm_value = (*l2norm)(upper)/(*l2norm)(lower);
	return norm_value*norm_value;
}

