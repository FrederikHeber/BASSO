/*
 * QuickAccessReferences.hpp
 *
 *  Created on: Dec 15, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_INVERSEPROBLEMS_QUICKACCESSREFERENCES_HPP_
#define MINIMIZATIONS_INVERSEPROBLEMS_QUICKACCESSREFERENCES_HPP_

#include "BassoConfig.h"

#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/types.hpp"

class LinearMapping;
class Mapping;
class NormedSpace;
class Norm;

/** This is a convenience class to contain quick and easy to remember
 * references to all instances of an inverse problem: Ax = y.
 *
 * Use like this:
 * \code
 * // instantiate with short name
 * QuickAccessReferences refs(_problem);
 * // and then use like this
 * SpaceElement_ptr_t xtemp = refs.SpaceX->createElement();
 * \endcode
 *
 */
struct QuickAccessReferences
{
	QuickAccessReferences(const InverseProblem_ptr_t &_problem) :
		SpaceX(*_problem->A->getSourceSpace()),
		DualSpaceX(*SpaceX.getDualSpace()),
		SpaceY(*_problem->A->getTargetSpace()),
		DualSpaceY(*SpaceY.getDualSpace()),
		NormX(*SpaceX.getNorm()),
		DualNormX(*DualSpaceX.getNorm()),
		NormY(*SpaceY.getNorm()),
		J_p(*SpaceX.getDualityMapping()),
		J_q(*DualSpaceX.getDualityMapping()),
		j_r(*SpaceY.getDualityMapping()),
		y(_problem->y),
		A(dynamic_cast<const LinearMapping &>(*_problem->A)),
		A_t(dynamic_cast<const LinearMapping &>(*A.getAdjointMapping()))
	{}

	// spaces
	//!> Source Space of linear Mapping \a A
	const NormedSpace & SpaceX;
	//!> Dual space to \a SpaceX
	const NormedSpace & DualSpaceX;
	//!> Target Space of linear Mapping \a A
	const NormedSpace & SpaceY;
	//!> Dual space to \a SpaceY
	const NormedSpace & DualSpaceY;

	// norms
	//!> Norm of \a SpaceX
	const Norm & NormX;
	//!> Norm of \a DualSpaceX
	const Norm & DualNormX;
	//!> Norm of \a SpaceY
	const Norm & NormY;

	// duality mappings
	//!> Duality Mapping from \a SpaceX to \a DualSpaceX
	const Mapping & J_p;
	//!> Duality Mapping from \a DualSpaceX to \a SpaceX
	const Mapping & J_q;
	//!> Duality Mapping from \a SpaceY to \a DualSpaceY
	const Mapping & j_r;

	// right-hand side
	//!> Right-hand side of the inverse problem
	const SpaceElement_ptr_t &y;

	// linear mappings
	//!> Linear mapping/operator to inverse problem
	const LinearMapping &A;
	//!> Adjoint mapping/operator of \a A
	const LinearMapping &A_t;
};


#endif /* MINIMIZATIONS_INVERSEPROBLEMS_QUICKACCESSREFERENCES_HPP_ */
