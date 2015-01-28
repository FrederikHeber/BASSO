/*
 * VectorSetter.hpp
 *
 *  Created on: Jan 28, 2015
 *      Author: heber
 */

#ifndef VECTORSETTER_HPP_
#define VECTORSETTER_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"

template <class vector_type>
struct VectorSetter
{
	static void set(
			const SpaceElement_ptr_t &_element,
			vector_type &_vector)
	{
		_vector = RepresentationAdvocate::get(_element);
	}
};

#endif /* VECTORSETTER_HPP_ */
