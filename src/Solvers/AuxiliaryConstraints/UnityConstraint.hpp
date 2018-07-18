/*
 * UnityConstraint.hpp
 *
 *  Created on: Nov 19, 2015
 *      Author: heber
 */

#ifndef SOLVERS_AUXILIARYCONSTRAINTS_UNITYCONSTRAINT_HPP_
#define SOLVERS_AUXILIARYCONSTRAINTS_UNITYCONSTRAINT_HPP_

#include "BassoConfig.h"

#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraint_impl.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"

struct UnityConstraint : public AuxiliaryConstraint_impl
{
	void operator()(SpaceElement_ptr_t &_vector) const
	{
		const double norm = _vector->Norm();
		*_vector *= 1./norm;
	}

	virtual ~UnityConstraint() {}
};




#endif /* SOLVERS_AUXILIARYCONSTRAINTS_UNITYCONSTRAINT_HPP_ */
