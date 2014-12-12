/*
 * Norm.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef NORM_HPP_
#define NORM_HPP_

#include "BassoConfig.h"

#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/VectorSpaceOperationCounts.hpp"

#include <Eigen/Dense>

/** This class defines the Norm interface.
 *
 */
class Norm
{
public:
	/** Constructor for class Norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 */
	Norm(const NormedSpace_ptr_t& _ref) :
		NormedSpaceRef(_ref)
	{}

	/** Virtual destructor.
	 *
	 */
	virtual ~Norm() {}

	/** Evaluates the norm for a given \a _element.
	 *
	 * This functions wraps the internal_operator() with an additional
	 * TimeKeeper object to measure how often the Norms are evaluated.
	 *
	 * @param _element element of the space, whose norm to evaluated
	 * @return norm of \a element
	 */
	const double operator()(const SpaceElement_ptr_t &_element) const
	{
		TIMEKEEPER(NormedSpaceRef->opcounts.VectorNorm);
		return internal_operator(_element);
	}

	/** Getter for the p value of a possible lp norm.
	 *
	 * @return p value: 0 - not an lp norm, else - p of lp norm
	 */
	virtual const double getPvalue() const
	{ return 0.; }

protected:
	//!> internal reference to the space this norm belongs to
	const NormedSpace_ptr_t NormedSpaceRef;

	/** Evaluates the norm for a given \a _element.
	 *
	 * @param _element element of the space, whose norm to evaluated
	 * @return norm of \a element
	 */
	virtual const double internal_operator(const SpaceElement_ptr_t &_element) const = 0;
};



#endif /* NORM_HPP_ */
