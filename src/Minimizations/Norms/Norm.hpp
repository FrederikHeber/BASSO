/*
 * Norm.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef NORM_HPP_
#define NORM_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

/** This class defines the Norm interface.
 *
 */
class Norm
{
public:
//	/** Constructor for class Norm.
//	 *
//	 * @param _ref reference to the space this norm is associated with
//	 */
//	Norm(const NormedSpace_ptr_t _ref) :
//		NormedSpaceRef(_ref)
//	{}

	/** Evaluates the norm for a given \a _element.
	 *
	 * @param _element element of the space, whose norm to evaluated
	 * @return norm of \a element
	 */
	virtual const double operator()(const SpaceElement_ptr_t &_element) const = 0;

//protected:
//	//!> internal reference to the space this norm belongs to
//	const NormedSpace_ptr_t NormedSpaceRef;
};



#endif /* NORM_HPP_ */
