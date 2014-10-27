/*
 * NormFactory.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormFactory.hpp"

#include <limits>

#include "Minimizations/Norms/L1Norm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/LInfinityNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"

Norm_ptr_t NormFactory::createInstance(const double _p)
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

