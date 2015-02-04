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

InverseProblem_ptr_t InverseProblemFactory::createLpInstance(
		const double _p,
		const double _powerx,
		const double _r,
		const double _powery,
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs)
{
	// first, create the spaces
	NormedSpace_ptr_t X =
			NormedSpaceFactory::createLpInstance(
					_matrix.outerSize(), _p, _powerx);
	NormedSpace_ptr_t Y =
			NormedSpaceFactory::createLpInstance(
					_matrix.innerSize(), _r, _powery);

	// then create the SpaceElement
	SpaceElement_ptr_t y = ElementCreator::create(Y, _rhs);

	// and the LinearMapping
	Mapping_ptr_t A = LinearMappingFactory::createInstance(X,Y,_matrix);

	InverseProblem_ptr_t instance(
			new InverseProblem(A,X,Y,y) );

	return instance;
}

InverseProblem_ptr_t InverseProblemFactory::createRegularizedL1Instance(
		const double _lambda,
		const double _powerx,
		const double _r,
		const double _powery,
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs)
{
	// first, create the spaces
	NormedSpace_ptr_t X =
			NormedSpaceFactory::createRegularizedL1Instance(
					_matrix.outerSize(), _lambda, _powerx);
	NormedSpace_ptr_t Y =
			NormedSpaceFactory::createLpInstance(
					_matrix.innerSize(), _r, _powery);

	// then create the SpaceElement
	SpaceElement_ptr_t y = ElementCreator::create(Y, _rhs);

	// and the LinearMapping
	Mapping_ptr_t A = LinearMappingFactory::createInstance(X,Y,_matrix);

	InverseProblem_ptr_t instance( new InverseProblem(A,X,Y,y) );

	return instance;
}
