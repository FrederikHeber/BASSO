/*
 * AuxiliaryConstraint_impl.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef SOLVERS_AUXILIARYCONSTRAINTS_AUXILIARYCONSTRAINT_IMPL_HPP_
#define SOLVERS_AUXILIARYCONSTRAINTS_AUXILIARYCONSTRAINT_IMPL_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

/** Interface definition for a stopping criterion.
 *
 * This defines only the function itself. In AuxiliaryConstraints we
 * wrap these "implementations" inside a hidden shared_ptr to allow
 * for combining these predicates.
 *
 */
struct AuxiliaryConstraint_impl
{
	/** Functor to check the respective stopping criterion.
	 *
	 * We have a broad interface, requesting all possible values
	 * where each of the stopping criteria needs only a few.
	 *
	 * @param _vector element to constrain
	 */
	virtual void operator()(SpaceElement_ptr_t &_vector) const = 0;
};

#endif /* SOLVERS_AUXILIARYCONSTRAINTS_AUXILIARYCONSTRAINT_IMPL_HPP_ */
