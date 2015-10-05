/*
 * SingularValueDecomposition.cpp
 *
 *  Created on: Oct 5, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include <cassert>

#include <Eigen/Dense>

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Mappings/SingularValueDecomposition.hpp"
#include "Minimizations/Mappings/SingularValueDecomposition_impl.hpp"
#include "Minimizations/Norms/Norm.hpp"

SpaceElement_ptr_t
SingularValueDecomposition::solve(SpaceElement_ptr_t _rhs) const
{
	// SVD only gives true solution for l2 norm
	assert( _rhs->getSpace()->getNorm()->getPvalue() == 2. );
	const Eigen::VectorXd &rhsvector =
			RepresentationAdvocate::get(_rhs);
	const Eigen::VectorXd solution = pimpl->svd.solve(rhsvector);
	return ElementCreator::create(_rhs->getSpace(), solution);
}
