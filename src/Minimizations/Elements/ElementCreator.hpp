/*
 * ElementCreator.hpp
 *
 *  Created on: Dec 11, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_ELEMENTS_ELEMENTCREATOR_HPP_
#define MINIMIZATIONS_ELEMENTS_ELEMENTCREATOR_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

/** This is just a helper class to hide the Eigen::... based
 * representation of SpaceElement's.
 */
struct ElementCreator
{
	static const SpaceElement_ptr_t create(
			const NormedSpace_ptr_t &_space,
			const Eigen::VectorXd &_vector);

	static const SpaceElement_ptr_t create(
			const NormedSpace &_space,
			const Eigen::VectorXd &_vector);
};


#endif /* MINIMIZATIONS_ELEMENTS_ELEMENTCREATOR_HPP_ */
