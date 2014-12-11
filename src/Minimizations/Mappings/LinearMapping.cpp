/*
 * LinearMapping.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "LinearMapping.hpp"

#include <cassert>

#include "Log/Logging.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

SpaceElement_ptr_t LinearMapping::operator()(
		const SpaceElement_ptr_t &_sourceelement
		) const
{
	assert( _sourceelement->getSpace() == SourceSpaceRef );
	SpaceElement_ptr_t targetelement = TargetSpaceRef->createElement();
	*targetelement =
			matrix * _sourceelement->getVectorRepresentation();
	return targetelement;
}

const Eigen::VectorXd LinearMapping::operator()(
		const Eigen::VectorXd &_sourcevector
		) const
{
	const Eigen::VectorXd targetvector = matrix * _sourcevector;
	return targetvector;
}

SpaceElement_ptr_t LinearMapping::operator*(const SpaceElement_ptr_t &_element) const
{
	assert( _element->getSpace() == SourceSpaceRef );
	SpaceElement_ptr_t newelement = TargetSpaceRef->createElement();
	*newelement =
			matrix * _element->getVectorRepresentation();
	return newelement;
}

const Eigen::VectorXd
LinearMapping::operator*(const Eigen::VectorXd &_vector) const
{
	const Eigen::VectorXd newvector = matrix * _vector;
	return newvector;
}

const Mapping_ptr_t LinearMapping::getAdjointMapping() const
{
	LinearMapping * adjoint = new LinearMapping(
				TargetSpaceRef->getDualSpace(),
				SourceSpaceRef->getDualSpace()
			);
	adjoint->matrix = matrix.transpose();
	return Mapping_ptr_t(adjoint);
}

const double LinearMapping::Norm() const
{
	BOOST_LOG_TRIVIAL(warning)
			<< "BEWARE: Is this calculating the right matrix norm?";
	return matrix.norm();
}

