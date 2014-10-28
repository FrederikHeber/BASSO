/*
 * NormedSpaceFactory.hpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */

#ifndef NORMEDSPACEFACTORY_HPP_
#define NORMEDSPACEFACTORY_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

/** This is a factory for spaces, creating all necessary member instances
 * and returning the NormedSpace instance.
 *
 */
struct NormedSpaceFactory
{
	/** Factory functor that creates the LpSpace instance.
	 *
	 * @param _dimension dimension of the space
	 * @param _p p value of the norm of the space
	 * @param _power power type of duality mapping's weight function
	 * @return NormedSpace instance according to parameters
	 */
	static NormedSpace_ptr_t createLpInstance(
			const unsigned int _dimension,
			const double _p,
			const double _power);

	/** Factory functor that creates the regularized L1 instance.
	 *
	 * @param _dimension dimension of the space
	 * @param _lambda regularization parameter
	 * @param _power power type of duality mapping's weight function
	 * @return NormedSpace instance according to parameters
	 */
	static NormedSpace_ptr_t createRegularizedL1Instance(
			const unsigned int _dimension,
			const double _lambda,
			const double _power);
};



#endif /* NORMEDSPACEFACTORY_HPP_ */
