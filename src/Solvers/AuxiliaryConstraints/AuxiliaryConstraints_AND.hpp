/*
 * AuxiliaryConstraints_AND.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef SOLVERS_AUXILIARYCONSTRAINTSS_AUXILIARYCONSTRAINT_AND_HPP_
#define SOLVERS_AUXILIARYCONSTRAINTSS_AUXILIARYCONSTRAINT_AND_HPP_

#include "BassoConfig.h"

#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints.hpp"

struct AuxiliaryConstraints_AND : public AuxiliaryConstraint_impl
{
	AuxiliaryConstraints_AND(
			const AuxiliaryConstraints::ptr_t &_left,
			const AuxiliaryConstraints::ptr_t &_right) :
				left(_left),
				right(_right)
	{}

	void operator()(SpaceElement_ptr_t &_vector) const
	{
		left->operator()(_vector);
		right->operator()(_vector);
	}

	virtual ~AuxiliaryConstraints_AND() {}

	const AuxiliaryConstraints::ptr_t left;
	const AuxiliaryConstraints::ptr_t right;
};



#endif /* SOLVERS_AUXILIARYCONSTRAINTSS_AUXILIARYCONSTRAINT_AND_HPP_ */
