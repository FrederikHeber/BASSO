/*
 * NormFactory.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormFactory.hpp"

#include <limits>

#include "Minimizations/Norms/IllegalNorm.hpp"
#include "Minimizations/Norms/L1Norm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/LInfinityNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/RegularizedL1Norm.hpp"

Norm_ptr_t NormFactory::createLpInstance(const double _p)
{
	Norm *NormY = NULL;
	if (_p == std::numeric_limits<double>::infinity()) {
		NormY = new LInfinityNorm;
	} else if (_p > 1.) {
		NormY = new LpNorm(_p);
	} else if (_p == 1.) {
		NormY = new L1Norm;
	}
	return Norm_ptr_t(NormY);
}

Norm_ptr_t NormFactory::createRegularizedL1Instance(
		const double _lambda)
{
	return Norm_ptr_t(new RegularizedL1Norm(_lambda));
}

Norm_ptr_t NormFactory::createIllegalInstance()
{
	return Norm_ptr_t(new IllegalNorm);
}
