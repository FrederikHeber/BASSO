/*
 * SearchspaceFactory.cpp
 *
 *  Created on: Nov 18, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SearchspaceFactory.hpp"

#include <cassert>

#include "Log/Logging.hpp"

#include "Minimizations/Minimizers/MinimizationExceptions.hpp"
#include "Minimizations/Minimizers/Searchspace/LastNSearchDirections.hpp"
#include "Minimizations/Minimizers/Searchspace/NemirovskyDirection.hpp"

// static entities
SearchspaceFactory::NameTypeMap_t SearchspaceFactory::NameTypeMap;
enum SearchspaceFactory::SearchspaceType
SearchspaceFactory::InstanceType = SearchspaceFactory::LastNDirections;

Searchspace::ptr_t SearchspaceFactory::create(
		const NormedSpace_ptr_t &_SearchDirectionSpace_ptr,
		const unsigned int _N,
		const OperationCounter<
			const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type,
			const Eigen::MatrixBase<Eigen::MatrixXd>&,
			const Eigen::MatrixBase<Eigen::VectorXd>&
			> &_MatrixVectorProduct_subspace,
		const OperationCounter<
			Eigen::internal::scalar_product_traits<typename Eigen::internal::traits<Eigen::VectorXd>::Scalar, typename Eigen::internal::traits<Eigen::VectorXd>::Scalar>::ReturnType,
			const Eigen::MatrixBase<Eigen::VectorXd>&,
			const Eigen::MatrixBase<Eigen::VectorXd>&
			> &_ScalarVectorProduct_subspace
		)
{
	Searchspace::ptr_t returninstance;
	switch (InstanceType) {
	case LastNDirections:
		returninstance.reset(new LastNSearchDirections(
				_SearchDirectionSpace_ptr,
				_N,
				_ScalarVectorProduct_subspace)
				);
		break;
	case Nemirovsky:
		if (_N != 2) {
			BOOST_LOG_TRIVIAL(error)
					<< "NemirovskyDirection always uses two search directions.";
			throw MinimizationIllegalValue_exception()
					<< MinimizationIllegalValue_name("N");
		}
		returninstance.reset(new NemirovskyDirection(
				_SearchDirectionSpace_ptr,
				_MatrixVectorProduct_subspace,
				_ScalarVectorProduct_subspace)
				);
		break;
	default:
		BOOST_LOG_TRIVIAL(error)
			<< "Current InstanceType is unknown.";
		assert(0);
	}
	return returninstance;
}

void
SearchspaceFactory::setCurrentType(const enum SearchspaceType _type)
{
	if (isValidType(_type))
		InstanceType = _type;
}

const std::string
SearchspaceFactory::getName(const enum SearchspaceType _type)
{
	fillNameTypeMap();
	if (isValidType(_type))
		return NameTypeMap.left.at(_type);
	else
		return NameTypeMap.left.at(InvalidType);
}

const enum SearchspaceFactory::SearchspaceType
SearchspaceFactory::getType(const std::string &_name)
{
	fillNameTypeMap();
	if (isValidName(_name))
		return NameTypeMap.right.at(_name);
	else
		return InvalidType;
}

bool
SearchspaceFactory::isValidName(const std::string &_name)
{
	fillNameTypeMap();
	NameTypeMap_t::right_const_iterator iter =
			NameTypeMap.right.find(_name);
	return ((iter != NameTypeMap.right.end())
			&& (iter->second != InvalidType));
}

bool
SearchspaceFactory::isValidType(const enum SearchspaceType _type)
{
	fillNameTypeMap();
	NameTypeMap_t::left_const_iterator iter =
			NameTypeMap.left.find(_type);
	return ((iter != NameTypeMap.left.end())
			&& (iter->first != InvalidType));
}

void
SearchspaceFactory::fillNameTypeMap()
{
	if (NameTypeMap.empty()) {
		NameTypeMap.insert(
				NameTypeMap_t::value_type( InvalidType, "") );
		NameTypeMap.insert(
				NameTypeMap_t::value_type( LastNDirections, "LastNDirections") );
		NameTypeMap.insert(
				NameTypeMap_t::value_type( Nemirovsky, "Nemirovsky") );
		// check that we inserted all
		assert( NameTypeMap.size() == (size_t)MAX_SearchspaceType );
	}
}
