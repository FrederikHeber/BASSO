/*
 * RelativeShrinkageCoefficient.hpp
 *
 *  Created on: Apr 25, 2017
 *      Author: heber
 */

#ifndef MINIMIZATIONS_NORMS_SPECIFICS_RELATIVESHRINKAGECOEFFICIENT_HPP_
#define MINIMIZATIONS_NORMS_SPECIFICS_RELATIVESHRINKAGECOEFFICIENT_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

/** This class calculates the relative shrinkage coefficient
 * that depends on the argument for the regularized l1-norm and
 * the RelativeShrinkageMapping.
 */
struct RelativeShrinkageCoefficient
{
	/** Calculates the relative shrinkage coefficient.
	 *
	 * @param _x
	 * @param _lambda regularization parameter
	 * @return coefficient alpha for the RelativeShrinkageMapping and RegularizedL1
	 */
	static const double
	get(
		const SpaceElement_ptr_t &_x,
		const double _lambda
			);
};



#endif /* MINIMIZATIONS_NORMS_SPECIFICS_RELATIVESHRINKAGECOEFFICIENT_HPP_ */
