/*
 * MappedNorm.hpp
 *
 *  Created on: Nov 9, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_NORMS_SPECIFICS_MAPPEDNORM_HPP_
#define MINIMIZATIONS_NORMS_SPECIFICS_MAPPEDNORM_HPP_

#include "BassoConfig.h"

#include <cassert>
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"

/** This class implements a linear mapped norm of the form
 * \f$ \tfrac 1 2 || Dx ||^2_{TV}\f$, where \f$ Dx \f$ is the
 * linear operator
 *
 */
class MappedNorm : public LpNorm
{
public:
	/** Constructor for class Norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _op linear mapping to modify argument before taking norm
	 */
	MappedNorm(
			const NormedSpace_weakptr_t& _ref,
			const LinearMapping &_op) :
		LpNorm(_ref, 2.),
		op(_op)
	{}

	bool isSmooth() const
	{ return true; }

protected:

	/** Evaluates the norm for a given \a _element.
	 *
	 * @param _element element of the space, whose norm to evaluated
	 * @return norm of \a element
	 */
	const double internal_operator(const SpaceElement_ptr_t &_x) const
	{
		assert( getSpace() == _x->getSpace() );
		return L2Norm::internal_operator(op(_x));
	}

private:
	//!> internal linear mapping to modify
	const LinearMapping op;
};


#endif /* MINIMIZATIONS_NORMS_SPECIFICS_MAPPEDNORM_HPP_ */
