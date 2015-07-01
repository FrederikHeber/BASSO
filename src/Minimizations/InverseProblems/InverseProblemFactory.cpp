/*
 * InverseProblemFactory.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InverseProblemFactory.hpp"

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

InverseProblem_ptr_t InverseProblemFactory::create(
		const std::string &_type_SpaceX,
		const args_t &_args_SpaceX,
		const std::string &_type_SpaceY,
		const args_t &_args_SpaceY,
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs)
{
	// first, create the spaces
	NormedSpace_ptr_t X = NormedSpaceFactory::create(
			_matrix.outerSize(), _type_SpaceX, _args_SpaceX);
	NormedSpace_ptr_t Y = NormedSpaceFactory::create(
			_matrix.innerSize(), _type_SpaceY, _args_SpaceY);

	// then create the SpaceElement
	SpaceElement_ptr_t y = ElementCreator::create(Y, _rhs);

	// and the LinearMapping
	Mapping_ptr_t A = LinearMappingFactory::createInstance(X,Y,_matrix);

	InverseProblem_ptr_t instance(
			new InverseProblem(A,X,Y,y) );

	return instance;
}

