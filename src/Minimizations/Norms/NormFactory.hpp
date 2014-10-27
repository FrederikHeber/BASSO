/*
 * NormFactory.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef NORMFACTORY_HPP_
#define NORMFACTORY_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

/** This factory instantiates the respective lp norm depending on the value
 * p.
 */
struct NormFactory
{
	/** Creates an lp norm according to the given \a _p.
	 *
	 * @param _p p value of the lp norm
	 * @return Norm instance
	 */
	static Norm_ptr_t createLpInstance(const double _p);

	/** Creates a regularized l1 norm.
	 *
	 * @param _lambda regularization parameter
	 * @return Norm instance
	 */
	static Norm_ptr_t createRegularizedL1Instance(
			const double _lambda);

	/** Creates an illegal norm.
	 *
	 * @return Norm instance
	 */
	static Norm_ptr_t createIllegalInstance();
};



#endif /* NORMFACTORY_HPP_ */
