/*
 * AuxiliaryConstraints.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: heber
 */

#ifndef SOLVERS_AUXILIARYCONSTRAINTS_AUXILIARYCONSTRAINTS_HPP_
#define SOLVERS_AUXILIARYCONSTRAINTS_AUXILIARYCONSTRAINTS_HPP_

#include "BassoConfig.h"

#include <boost/shared_ptr.hpp>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Solvers/AuxiliaryConstraints/AuxiliaryConstraint_impl.hpp"

/** Interface definition for a stopping criterion.
 *
 */
struct AuxiliaryConstraints
{
	//!> typedef for shared ptr containing AuxiliaryConstraints
	typedef boost::shared_ptr<AuxiliaryConstraint_impl> ptr_t;

	/** Functor to check the respective stopping criterion.
	 *
	 * We have a broad interface, requesting all possible values
	 * where each of the stopping criteria needs only a few.
	 *
	 * @param _vector element to constrain
	 */
	void operator()(
			SpaceElement_ptr_t &_vector) const;

private:
	//!> internal predicate evaluating the stopping criterion or a combination
	const ptr_t internal_impl;
};

// boolean operators to combine stopping criteria
AuxiliaryConstraints::ptr_t operator&&(
		const AuxiliaryConstraints::ptr_t &_a,
		const AuxiliaryConstraints::ptr_t &_b);

#endif /* SOLVERS_AUXILIARYCONSTRAINTS_AUXILIARYCONSTRAINTS_HPP_ */
