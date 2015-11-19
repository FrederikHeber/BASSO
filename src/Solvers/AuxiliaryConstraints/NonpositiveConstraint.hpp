/*
 * NonpositiveConstraint.hpp
 *
 *  Created on: Nov 19, 2015
 *      Author: heber
 */

#ifndef SOLVERS_AUXILIARYCONSTRAINTS_NONPOSITIVECONSTRAINT_HPP_
#define SOLVERS_AUXILIARYCONSTRAINTS_NONPOSITIVECONSTRAINT_HPP_

#include "BassoConfig.h"

#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraint_impl.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"

struct NonpositiveConstraint : public AuxiliaryConstraint_impl
{
	void operator()(SpaceElement_ptr_t &_vector) const
	{
		const SpaceElement_ptr_t tempvector = (-1.*_vector)->getAbsVector();
		*_vector = .5*(_vector - tempvector);
	}
};




#endif /* SOLVERS_AUXILIARYCONSTRAINTS_NONPOSITIVECONSTRAINT_HPP_ */
