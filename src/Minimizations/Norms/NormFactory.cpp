/*
 * NormFactory.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormFactory.hpp"

#include <limits>

#include "Minimizations/Norms/DualRegularizedL1Norm.hpp"
#include "Minimizations/Norms/IllegalNorm.hpp"
#include "Minimizations/Norms/L1Norm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/LInfinityNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/RegularizedL1Norm.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

Norm_ptr_t NormFactory::createLpInstance(
		const double _p)
{
	Norm *NormY = NULL;
	if (_p == std::numeric_limits<double>::infinity()) {
		NormY = new LInfinityNorm(NormedSpaceFactory::DummySpace);
	} else if (_p > 1.) {
		NormY = new LpNorm(NormedSpaceFactory::DummySpace, _p);
	} else if (_p == 1.) {
		NormY = new L1Norm(NormedSpaceFactory::DummySpace);
	} else
		throw NormIllegalValue_exception()
			<< NormIllegalValue_name("p");

	return Norm_ptr_t(NormY);
}

Norm_ptr_t NormFactory::createLpInstance(
		const NormedSpace_ptr_t& _ref,
		const double _p)
{
	Norm *NormY = NULL;
	if (_p == std::numeric_limits<double>::infinity()) {
		NormY = new LInfinityNorm(_ref);
	} else if (_p > 1.) {
		NormY = new LpNorm(_ref, _p);
	} else if (_p == 1.) {
		NormY = new L1Norm(_ref);
	}
	return Norm_ptr_t(NormY);
}

Norm_ptr_t NormFactory::createRegularizedL1Instance(
		const NormedSpace_ptr_t& _ref,
		const double _lambda)
{
	return Norm_ptr_t(new RegularizedL1Norm(_ref, _lambda));
}

Norm_ptr_t NormFactory::createDualRegularizedL1Instance(
		const NormedSpace_ptr_t& _ref,
		const double _lambda)
{
	return Norm_ptr_t(new DualRegularizedL1Norm(_ref, _lambda));
}

Norm_ptr_t NormFactory::createIllegalInstance(
		const NormedSpace_ptr_t& _ref)
{
	return Norm_ptr_t(new IllegalNorm(_ref));
}
