/*
 * IllegalDualityMapping.cpp
 *
 *  Created on: Oct 14, 2014
 *      Author: heber
 */

#include "IllegalDualityMapping.hpp"

#include "Minimizations/MinimizationExceptions.hpp"

IllegalDualityMapping::IllegalDualityMapping() :
	LpDualityMapping(1.)
{}

const Eigen::VectorXd
IllegalDualityMapping::operator()(
		const Eigen::VectorXd &_x,
		const double _power) const
{
	// we just throw as this function must now be called
	throw MinimizationIllegalValue_exception()
		<< MinimizationIllegalValue_name("illegally called");
}


