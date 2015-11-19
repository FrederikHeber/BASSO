/*
 * AuxiliaryConstraints.cpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "AuxiliaryConstraints.hpp"

#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraints_AND.hpp"

void AuxiliaryConstraints::operator()(
		SpaceElement_ptr_t &_vector) const
{
	internal_impl->operator()(_vector);
}

AuxiliaryConstraints::ptr_t operator&&(
		const AuxiliaryConstraints::ptr_t &_a,
		const AuxiliaryConstraints::ptr_t &_b)
{
	AuxiliaryConstraints::ptr_t criterion(new AuxiliaryConstraints_AND(_a, _b));
	return criterion;
}

