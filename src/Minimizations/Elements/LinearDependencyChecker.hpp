/*
 * LinearDependencyChecker.hpp
 *
 *  Created on: Feb 10, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_ELEMENTS_LINEARDEPENDENCYCHECKER_HPP_
#define MINIMIZATIONS_ELEMENTS_LINEARDEPENDENCYCHECKER_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

#include <vector>

/** This functor checks for a set of SpaceElements whether the representations
 * are linear dependent or not.
 */
struct LinearDependencyChecker
{
	//!> typedef for the set of SpaceElements
	typedef std::vector<SpaceElement_ptr_t> vectors_t;

	/** Checks whether in the given set there are two vectors that are
	 * linear dependent.
	 *
	 * \note We do this by constructing a matrix and performing a Gaussian
	 * Elimination on it.
	 *
	 * @param _vectors set of SpaceElements to check
	 * @return true - at least one pair is linear dependent, false - else
	 */
	bool operator()(const vectors_t &_vectors);
};


#endif /* MINIMIZATIONS_ELEMENTS_LINEARDEPENDENCYCHECKER_HPP_ */
