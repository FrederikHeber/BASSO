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
	static Norm_ptr_t createInstance(const double _p);
};



#endif /* NORMFACTORY_HPP_ */
