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

class L1NormUnitTest;
class LInfinityNormUnitTest;
class LpNormUnitTest;

/** This factory instantiates the respective lp norm depending on the value
 * p.
 */
class NormFactory
{
	//!> grant L1Norm unit test access to "Space-less" creator function
	friend class L1NormUnitTest;
	//!> grant LInfinityNorm unit test access to "Space-less" creator function
	friend class LInfinityNormUnitTest;
	//!> grant LpNorm unit test access to "Space-less" creator function
	friend class LpNormUnitTest;

	/** Creates an lp norm according to the given \a _p.
	 *
	 * @param _p p value of the lp norm
	 * @return Norm instance
	 */
	static Norm_ptr_t createLpInstance(const double _p);

public:
	/** Creates an lp norm according to the given \a _p.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _p p value of the lp norm
	 * @return Norm instance
	 */
	static Norm_ptr_t createLpInstance(
			const NormedSpace_weakptr_t _ref,
			const double _p);

	/** Creates a regularized l1 norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _lambda regularization parameter
	 * @return Norm instance
	 */
	static Norm_ptr_t createRegularizedL1Instance(
			const NormedSpace_weakptr_t _ref,
			const double _lambda);

	/** Creates a dual norm to regularized l1 norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @param _lambda regularization parameter
	 * @return Norm instance
	 */
	static Norm_ptr_t createDualRegularizedL1Instance(
			const NormedSpace_weakptr_t _ref,
			const double _lambda);

	/** Creates an illegal norm.
	 *
	 * @param _ref reference to the space this norm is associated with
	 * @return Norm instance
	 */
	static Norm_ptr_t createIllegalInstance(
			const NormedSpace_weakptr_t _ref);
};



#endif /* NORMFACTORY_HPP_ */
