/*
 * InverseProblemFactory.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "InverseProblemFactory.hpp"

#include <boost/assign.hpp>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"

using namespace boost::assign;

InverseProblem_ptr_t InverseProblemFactory::createLpInstance(
		const double _p,
		const double _powerx,
		const double _r,
		const double _powery,
		const Eigen::MatrixXd &_matrix,
		const Eigen::VectorXd &_rhs)
{
	// first, create the spaces
	NormedSpace_ptr_t X;
	{
		NormedSpaceFactory::args_t args;
		args += boost::any(_p), boost::any(_powerx);
		X = NormedSpaceFactory::create(
						_matrix.outerSize(), "lp", args);
	}
	NormedSpace_ptr_t Y;
	{
		NormedSpaceFactory::args_t args;
		args += boost::any(_r), boost::any(_powery);
		Y = NormedSpaceFactory::create(
						_matrix.innerSize(), "lp", args);
	}

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
	NormedSpace_ptr_t X;
	{
		NormedSpaceFactory::args_t args;
		args += boost::any(_lambda), boost::any(_powerx);
		X = NormedSpaceFactory::create(
						_matrix.outerSize(), "regularized_l1", args);
	}
	NormedSpace_ptr_t Y;
	{
		NormedSpaceFactory::args_t args;
		args += boost::any(_r), boost::any(_powery);
		Y = NormedSpaceFactory::create(
						_matrix.innerSize(), "lp", args);
	}

	// then create the SpaceElement
	SpaceElement_ptr_t y = ElementCreator::create(Y, _rhs);

	// and the LinearMapping
	Mapping_ptr_t A = LinearMappingFactory::createInstance(X,Y,_matrix);

	InverseProblem_ptr_t instance( new InverseProblem(A,X,Y,y) );

	return instance;
}
