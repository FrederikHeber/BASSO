/*
 * SoftThresholdingL1Norm.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef SOFTTHRESHOLDINGL1NORM_HPP_
#define SOFTTHRESHOLDINGL1NORM_HPP_

#include "BassoConfig.h"

#include <cassert>
#include "Minimizations/Norms/Specifics/RegularizedL1Norm.hpp"

/** This class is just a renaming of the regularized l1 norm such that
 * we have a unique dual norm to it.
 *
 */
struct SoftThresholdingL1Norm : public RegularizedL1Norm
{
public:
	/** Constructor for class Norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _lambda regularization parameter
	 */
	SoftThresholdingL1Norm(
			const NormedSpace_weakptr_t& _ref,
			const double _lambda = 0.1) :
		RegularizedL1Norm(_ref, _lambda)
	{}
};

#endif /* SOFTTHRESHOLDINGL1NORM_HPP_ */
