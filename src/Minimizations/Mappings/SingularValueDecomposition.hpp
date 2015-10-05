/*
 * SingularValueDecomposition.hpp
 *
 *  Created on: Oct 5, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MAPPINGS_SINGULARVALUEDECOMPOSITION_HPP_
#define MINIMIZATIONS_MAPPINGS_SINGULARVALUEDECOMPOSITION_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

struct SingularValueDecomposition
{
	SingularValueDecomposition(
			SingularValueDecomposition_impl_ptr_t &_pimpl) :
		pimpl(_pimpl)
	{}

	SpaceElement_ptr_t solve(SpaceElement_ptr_t _rhs) const;

private:
	SingularValueDecomposition_impl_ptr_t &pimpl;
};



#endif /* MINIMIZATIONS_MAPPINGS_SINGULARVALUEDECOMPOSITION_HPP_ */
