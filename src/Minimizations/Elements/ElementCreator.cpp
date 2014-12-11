/*
 * ElementCreator.cpp
 *
 *  Created on: Dec 11, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "ElementCreator.hpp"

const SpaceElement_ptr_t ElementCreator::create(
		const NormedSpace_ptr_t &_space,
		const Eigen::VectorXd &_vector)
{
	return create(*_space, _vector);
}

const SpaceElement_ptr_t ElementCreator::create(
		const NormedSpace &_space,
		const Eigen::VectorXd &_vector)
{
	SpaceElement_ptr_t returnelement = _space.createElement();
	*returnelement = _vector;
	return returnelement;
}
