/*
 * DualityMappingFactory.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */


#include <Minimizations/Mappings/DualityMappingFactory.hpp>
#include "BassoConfig.h"

#include <cassert>
#include <limits>

#include <boost/bind.hpp>

#include "Log/Logging.hpp"
#include "Math/Helpers.hpp"
#include "Minimizations/Mappings/IllegalDualityMapping.hpp"
#include "Minimizations/Mappings/MappingExceptions.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/L1DualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LpDualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/LInfinityDualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/RegularizedL1DualityMapping.hpp"
#include "Minimizations/Mappings/Specifics/RelativeShrinkageMapping.hpp"
#include "Minimizations/Mappings/Specifics/SoftThresholdingMapping.hpp"

const DualityMappingFactory::TokenCreatorMap_t DualityMappingFactory::getMap()
{
	TokenCreatorMap_t TokenCreatorMap;
	// construct TokenCreatorMap_t on constructing singleton instance
	TokenCreatorMap["lp"] =
			boost::bind(&DualityMappingFactory::createPowerTypeInstance, _1, _2);
	TokenCreatorMap["dual_lp"] =
			boost::bind(&DualityMappingFactory::createDualPowerTypeInstance, _1, _2);
	// regularized_1 is not smooth, we use reularized l1 mapping
	TokenCreatorMap["regularized_l1"] =
			boost::bind(&DualityMappingFactory::createRegularizedL1Instance, _1, _2);
	// dual of regularized_1 is smooth, is relative shrinkage
	TokenCreatorMap["dual_regularized_l1"] =
			boost::bind(&DualityMappingFactory::createRelativeShrinkrageInstance, _1, _2);
	// regularized_1 is not smooth, hence illegal mapping
	TokenCreatorMap["softthresholding_l1"] =
			boost::bind(&DualityMappingFactory::createRegularizedL1Instance, _1, _2);
	// dual of regularized_1 is smooth, is relative shrinkages
	TokenCreatorMap["dual_softthresholding_l1"] =
			boost::bind(&DualityMappingFactory::createSoftThresholdingInstance, _1, _2);
	TokenCreatorMap["relativeshrinkage_l1"] =
			boost::bind(&DualityMappingFactory::createRegularizedL1Instance, _1, _2);
	// dual of regularized_1 is smooth, is relative shrinkage
	TokenCreatorMap["dual_relativeshrinkage_l1"] =
			boost::bind(&DualityMappingFactory::createRelativeShrinkrageInstance, _1, _2);
	return TokenCreatorMap;
}

Mapping_ptr_t DualityMappingFactory::create(
		const std::string &_token,
		const NormedSpace_weakptr_t _space,
		const args_t &_args)
{
	assert( isValidType(_token) );
	const TokenCreatorMap_t creatormap = DualityMappingFactory::getMap();
	TokenCreatorMap_t::const_iterator iter =
			creatormap.find(_token);
	return (iter->second)(_space, _args);
}

bool DualityMappingFactory::isValidType(
		const std::string &_token
		)
{
	const TokenCreatorMap_t creatormap = DualityMappingFactory::getMap();
	TokenCreatorMap_t::const_iterator iter =
			creatormap.find(_token);
	return (iter != creatormap.end());
}


template <typename T>
const double getNthArgumentAs(
		const DualityMappingFactory::args_t &_args,
		const size_t _N)
{
	if(_args.size() < _N)
		throw MappingIllegalValue_exception()
			<< MappingIllegalValue_name("Not enough arguments given.");
    try {
        const T value = boost::any_cast<const T>(_args[_N]);
//    	LOG(info, "getNthArgumentAs of args " << "[" << _N << "] is " << value;
    	return value;
    } catch(const boost::bad_any_cast &) {
		throw MappingIllegalValue_exception()
			<< MappingIllegalValue_name("Nth argument is not of desired type");
    }
}

static PowerTypeDualityMapping* createPowerTypeDualityMapping(
		const NormedSpace_weakptr_t _NormedSpaceRef,
		const double _power
		)
{
	PowerTypeDualityMapping *mapping = NULL;
	const double p = NormedSpace_ptr_t(_NormedSpaceRef)->getNorm()->getPvalue();
	if (p == std::numeric_limits<double>::infinity()) {
		mapping = new LInfinityDualityMapping(_NormedSpaceRef, _power);
	} else if (p > 1.) {
		mapping = new LpDualityMapping(_NormedSpaceRef, _power);
	} else if (p == 1.) {
		mapping = new L1DualityMapping(_NormedSpaceRef, _power);
	}
	return mapping;
}

Mapping_ptr_t
DualityMappingFactory::createPowerTypeInstance(
		const NormedSpace_weakptr_t _NormedSpaceRef,
		const args_t &_args)
{
	const double power = getNthArgumentAs<double>(_args, 0);
	return Mapping_ptr_t(
			createPowerTypeDualityMapping(_NormedSpaceRef, power));
}

Mapping_ptr_t
DualityMappingFactory::createDualPowerTypeInstance(
		const NormedSpace_weakptr_t _NormedSpaceRef,
		const args_t &_args)
{
	const double power = getNthArgumentAs<double>(_args, 0);
	const double qpower = Helpers::ConjugateValue(power);
	return Mapping_ptr_t(
			createPowerTypeDualityMapping(_NormedSpaceRef, qpower));
}

Mapping_ptr_t
DualityMappingFactory::createIllegalInstance(
		const NormedSpace_weakptr_t _NormedSpaceRef,
		const args_t &_args)
{
	DualityMapping *mapping = new IllegalDualityMapping;
	return Mapping_ptr_t(mapping);
}

Mapping_ptr_t
DualityMappingFactory::createSoftThresholdingInstance(
		const NormedSpace_weakptr_t _NormedSpaceRef,
		const args_t &_args)
{
	const double lambda = getNthArgumentAs<double>(_args, 0);
	DualityMapping *mapping = new SoftThresholdingMapping(
			_NormedSpaceRef, lambda);
	return Mapping_ptr_t(mapping);
}

Mapping_ptr_t
DualityMappingFactory::createRelativeShrinkrageInstance(
		const NormedSpace_weakptr_t _NormedSpaceRef,
		const args_t &_args)
{
	const double lambda = getNthArgumentAs<double>(_args, 0);
	DualityMapping *mapping = new RelativeShrinkageMapping(
			_NormedSpaceRef, lambda);
	return Mapping_ptr_t(mapping);
}

Mapping_ptr_t
DualityMappingFactory::createRegularizedL1Instance(
		const NormedSpace_weakptr_t _NormedSpaceRef,
		const args_t &_args)
{
	const double lambda = getNthArgumentAs<double>(_args, 0);
	DualityMapping *mapping = new RegularizedL1DualityMapping(
			_NormedSpaceRef, lambda);
	return Mapping_ptr_t(mapping);
}
