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

struct VectorSetter
{
	/** Setter for vector and matrix cols/rows.
	 *
	 * Use like this be first defining an Eigen::Ref to your object:
	 * \begincode
	 * const Eigen::Ref<Eigen::MatrixXd::ColXpr> col = matrix.col(i);
	 * VectorSetter::set(col, ...);
	 * \endcode
	 *
	 * @param _vector
	 * @param _element
	 */
	template <class T>
	static void set(
			Eigen::Ref<T> _vector,
			const SpaceElement_ptr_t &_element)
	{
		_vector = RepresentationAdvocate::get(_element);
	}
};

#endif /* VECTORSETTER_HPP_ */
