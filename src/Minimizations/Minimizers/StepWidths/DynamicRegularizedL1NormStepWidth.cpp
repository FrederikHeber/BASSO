/*
 * DynamicRegularizedL1NormStepWidth.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "DynamicRegularizedL1NormStepWidth.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/NormFactory.hpp"

DynamicRegularizedL1NormStepWidth::DynamicRegularizedL1NormStepWidth(
		const InverseProblem_ptr_t &_problem,
		const Mapping_ptr_t &_J_r) :
	problem(_problem),
	J_r(_J_r),
	A_adjoint(_problem->A->getAdjointMapping()),
	l2norm_Y(NormFactory::create(
			"lp",
			problem->A->getTargetSpace(),
			NormFactory::args_t(1, boost::any(2.)))),
	l2norm_DualX(NormFactory::create(
			"dual_lp",
			problem->x->getSpace()->getDualSpace(),
			NormFactory::args_t(1, boost::any(2.))))
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
	const LinearMapping &A = static_cast<const LinearMapping &>(*problem->A);
	const LinearMapping &A_t = static_cast<const LinearMapping &>(*A_adjoint);
	const SpaceElement_ptr_t upper = A * _solution - problem->y;
	const SpaceElement_ptr_t lower = A_t * (*J_r)(upper);
	const double norm_value = (*l2norm_Y)(upper)/(*l2norm_DualX)(lower);
	return norm_value*norm_value;
}

