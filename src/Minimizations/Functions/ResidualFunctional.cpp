/*
 * ResidualFunctional.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "ResidualFunctional.hpp"

double ResidualFunctional::operator()(double _arg) const
{
	SpaceElement_ptr_t dual_solution = dualx->getSpace()->createElement();
	*dual_solution = dualx;
	*dual_solution -= _arg * u;
	*problem->x =
			(*problem->x->getSpace()->getDualSpace()->getDualityMapping())(
					dual_solution);
	// calculate residual at candidate position (goes into hx[0])
	return mininizer.calculateResidual( problem, residualvector );
}

