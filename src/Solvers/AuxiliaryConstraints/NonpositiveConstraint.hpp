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

#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"

struct NonpositiveConstraint : public AuxiliaryConstraint_impl
{
	void operator()(SpaceElement_ptr_t &_vector) const
	{
		RepresentationAdvocate::set(
				_vector,
				.5*((-1.*RepresentationAdvocate::get(_vector)).array().abs()-
						RepresentationAdvocate::get(_vector).array()));
	}
};




#endif /* SOLVERS_AUXILIARYCONSTRAINTS_NONPOSITIVECONSTRAINT_HPP_ */
