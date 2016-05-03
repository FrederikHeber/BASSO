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
	*dual_solution = dualx;
	*dual_solution -= _arg * u;
	*problem->x =
			(*problem->SourceSpace->getDualSpace()->getDualityMapping())(
					dual_solution);
	// calculate residual at candidate position (goes into hx[0])
	return residualizer( problem, residualvector );
}

