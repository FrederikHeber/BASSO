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
	const Eigen::VectorXd upper =
			dynamic_cast<const LinearMapping &>(*problem->A).getMatrixRepresentation()
			* _solution->getVectorRepresentation()
			- problem->y->getVectorRepresentation();
	Eigen::VectorXd lower =
			dynamic_cast<const LinearMapping &>(*problem->A).getMatrixRepresentation()
			* _solution->getVectorRepresentation()
			- problem->y->getVectorRepresentation();
	lower = dynamic_cast<const LinearMapping &>(*A_adjoint).getMatrixRepresentation()
			* lower;
	return upper.squaredNorm()/lower.squaredNorm();
}

