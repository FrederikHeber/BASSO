/*
 * InverseProblemFactory.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InverseProblemFactory.hpp"

#include <boost/bind.hpp>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

InverseProblem_ptr_t InverseProblemFactory::create_internal(
		const std::string &_type_SpaceX,
		const args_t &_args_SpaceX,
		const std::string &_type_SpaceY,
		const args_t &_args_SpaceY,
		const Eigen::VectorXd &_rhs,
		const unsigned int _outerSize,
		const unsigned int _innerSize,
		const MappingCreator_t &_creator
		)
{
	// first, create the spaces
	NormedSpace_ptr_t X = NormedSpaceFactory::create(
			_outerSize, _type_SpaceX, _args_SpaceX);
	NormedSpace_ptr_t Y = NormedSpaceFactory::create(
			_innerSize, _type_SpaceY, _args_SpaceY);

	// then create the SpaceElement
	SpaceElement_ptr_t y = ElementCreator::create(Y, _rhs);

	// and the LinearMapping
	Mapping_ptr_t A = _creator(X,Y);

	InverseProblem_ptr_t instance(
			new InverseProblem(A,X,Y,y) );

	return instance;
}

InverseProblem_ptr_t InverseProblemFactory::create(
		const std::string &_type_SpaceX,
		const args_t &_args_SpaceX,
		const std::string &_type_SpaceY,
		const args_t &_args_SpaceY,
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs)
{
	const MappingCreator_t creator = boost::bind(&LinearMappingFactory::createInstance,
			_1, _2, boost::cref(_matrix), false);
	return create_internal(
			_type_SpaceX, _args_SpaceX,
			_type_SpaceY, _args_SpaceY,
			_rhs,
			_matrix.outerSize(), _matrix.innerSize(),
			creator);
}

InverseProblem_ptr_t InverseProblemFactory::createFromTwoFactors(
		const std::string &_type_SpaceX,
		const args_t &_args_SpaceX,
		const std::string &_type_SpaceY,
		const args_t &_args_SpaceY,
		const Eigen::MatrixXd &_matrix_first_factor,
		const Eigen::MatrixXd &_matrix_second_factor,
		const Eigen::VectorXd &_rhs)
{
	const MappingCreator_t creator = boost::bind(&LinearMappingFactory::createTwoFactorInstance,
			_1, _2, boost::cref(_matrix_first_factor), boost::cref(_matrix_second_factor), false);
	return create_internal(
			_type_SpaceX, _args_SpaceX,
			_type_SpaceY, _args_SpaceY,
			_rhs,
			_matrix_first_factor.innerSize(), _matrix_second_factor.outerSize(),
			creator);
}
