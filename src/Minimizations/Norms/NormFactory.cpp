/*
 * NormFactory.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NormFactory.hpp"

#include <limits>

#include <boost/bind.hpp>

#include "Log/Logging.hpp"
#include "Math/Helpers.hpp"
#include "Minimizations/Norms/IllegalNorm.hpp"
#include "Minimizations/Norms/L1Norm.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/LInfinityNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Norms/Specifics/DualRelativeShrinkageL1Norm.hpp"
#include "Minimizations/Norms/Specifics/DualSoftThresholdingL1Norm.hpp"
#include "Minimizations/Norms/Specifics/RelativeShrinkageL1Norm.hpp"
#include "Minimizations/Norms/Specifics/SoftThresholdingL1Norm.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

const NormFactory::TokenCreatorMap_t NormFactory::getMap()
{
	TokenCreatorMap_t TokenCreatorMap;
	// construct TokenCreatorMap_t on constructing singleton instance
//	TokenCreatorMap["illegal"] =
//			boost::bind(&NormFactory::createIllegalInstance, _1, _2);
	TokenCreatorMap["lp"] =
			boost::bind(&NormFactory::createLpInstance, _1, _2);
	TokenCreatorMap["dual_lp"] =
			boost::bind(&NormFactory::createDualLpInstance, _1, _2);
	TokenCreatorMap["regularized_l1"] =
			boost::bind(&NormFactory::createRelativeShrinkageL1Instance, _1, _2);
	TokenCreatorMap["dual_regularized_l1"] =
			boost::bind(&NormFactory::createDualRelativeShrinkageL1Instance, _1, _2);
	TokenCreatorMap["relativeshrinkage_l1"] =
			boost::bind(&NormFactory::createRelativeShrinkageL1Instance, _1, _2);
	TokenCreatorMap["dual_relativeshrinkage_l1"] =
			boost::bind(&NormFactory::createDualRelativeShrinkageL1Instance, _1, _2);
	TokenCreatorMap["softthresholding_l1"] =
			boost::bind(&NormFactory::createSoftThresholdingL1Instance, _1, _2);
	TokenCreatorMap["dual_softthresholding_l1"] =
			boost::bind(&NormFactory::createDualSoftThresholdingL1Instance, _1, _2);
	return TokenCreatorMap;
}

Norm_ptr_t NormFactory::create(
		const std::string &_token,
		const NormedSpace_weakptr_t _space,
		const args_t &_args)
{
	const TokenCreatorMap_t creatormap = NormFactory::getMap();
	TokenCreatorMap_t::const_iterator iter =
			creatormap.find(_token);
	assert( iter != creatormap.end() );
	return (iter->second)(_space, _args);
}

bool NormFactory::isValidType(
		const std::string &_token
		)
{
	const TokenCreatorMap_t creatormap = NormFactory::getMap();
	TokenCreatorMap_t::const_iterator iter =
			creatormap.find(_token);
	return (iter != creatormap.end());
}

template <typename T>
const double getNthArgumentAs(
		const NormFactory::args_t &_args,
		const size_t _N)
{
	if(_args.size() < _N)
		throw NormIllegalValue_exception()
			<< NormIllegalValue_name("Not enough arguments given.");
    try {
        const T value = boost::any_cast<const T>(_args[_N]);
//    	LOG(info, "getNthArgumentAs of args " << "[" << _N << "] is " << value);
    	return value;
    } catch(const boost::bad_any_cast &) {
		throw NormIllegalValue_exception()
			<< NormIllegalValue_name("Nth argument is not of desired type");
    }
}

static Norm_ptr_t createLpNorm(
		const NormedSpace_weakptr_t _ref,
		const double p)
{
	Norm *NormX = NULL;
	if (p == std::numeric_limits<double>::infinity()) {
		NormX = new LInfinityNorm(_ref);
	} else if (p > 1.) {
		NormX = new LpNorm(_ref, p);
	} else if (p == 1.) {
		NormX = new L1Norm(_ref);
	} else
		throw NormIllegalValue_exception()
			<< NormIllegalValue_name("p must be in {0, [1,\\infty]}");

	return Norm_ptr_t(NormX);
}

Norm_ptr_t NormFactory::createLpInstance(
		const double _p)
{
	// give warning for very small or large lp norms
	if ((_p < 1.01) || (_p > 100)) {
		LOG(warning, "The specified p value of " << _p
			<< " is very small and will most likely lead to numerical instabilities." << "It is suggested to use a regularized l1 or linf norm.");
	}

	return createLpNorm(NormedSpaceFactory::getDummySpace(), _p);
}

Norm_ptr_t NormFactory::createLpInstance(
		const NormedSpace_weakptr_t _ref,
		const args_t &_args)
{
	const double p = getNthArgumentAs<double>(_args, 0);
	return createLpNorm(_ref, p);
}

Norm_ptr_t NormFactory::createDualLpInstance(
		const NormedSpace_weakptr_t _ref,
		const args_t &_args)
{
	const double p = getNthArgumentAs<double>(_args, 0);
	const double q = Helpers::ConjugateValue(p);
	return createLpNorm(_ref, q);
}

Norm_ptr_t NormFactory::createRelativeShrinkageL1Instance(
		const NormedSpace_weakptr_t _ref,
		const args_t &_args)
{
	const double lambda = getNthArgumentAs<double>(_args, 0);
	return Norm_ptr_t(new RelativeShrinkageL1Norm(_ref, lambda));
}

Norm_ptr_t NormFactory::createDualRelativeShrinkageL1Instance(
		const NormedSpace_weakptr_t _ref,
		const args_t &_args)
{
	const double lambda = getNthArgumentAs<double>(_args, 0);
	return Norm_ptr_t(new DualRelativeShrinkageL1Norm(_ref, lambda));
}

Norm_ptr_t NormFactory::createSoftThresholdingL1Instance(
		const NormedSpace_weakptr_t _ref,
		const args_t &_args)
{
	const double lambda = getNthArgumentAs<double>(_args, 0);
	return Norm_ptr_t(new SoftThresholdingL1Norm(_ref, lambda));
}

Norm_ptr_t NormFactory::createDualSoftThresholdingL1Instance(
		const NormedSpace_weakptr_t _ref,
		const args_t &_args)
{
	const double lambda = getNthArgumentAs<double>(_args, 0);
	return Norm_ptr_t(new DualSoftThresholdingL1Norm(_ref, lambda));
}

Norm_ptr_t NormFactory::createIllegalInstance(
		const NormedSpace_weakptr_t _ref,
		const args_t &_args)
{
	return Norm_ptr_t(new IllegalNorm(_ref));
}
