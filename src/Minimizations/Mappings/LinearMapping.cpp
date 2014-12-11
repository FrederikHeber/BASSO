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

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

const SpaceElement_ptr_t LinearMapping::operator()(
		const SpaceElement_ptr_t &_sourceelement
		) const
{
	assert( _sourceelement->getSpace() == SourceSpaceRef );
	SpaceElement_ptr_t targetelement =
			ElementCreator::create(
					TargetSpaceRef,
					matrix * RepresentationAdvocate::get(_sourceelement));
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
	SpaceElement_ptr_t targetelement =
			ElementCreator::create(
					TargetSpaceRef,
					matrix * RepresentationAdvocate::get(_element));
	return targetelement;
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

const double LinearMapping::MutualCoherence() const
{
	double mutual_coherence = 0.;
	for (unsigned int i=0;i<getSourceSpace()->getDimension();++i) {
		for (unsigned int j=i+1;j<getSourceSpace()->getDimension();++j) {
			const Eigen::VectorXd col_i = matrix.col(i);
			const Eigen::VectorXd col_j = matrix.col(j);
			const double col_i_norm = col_i.norm();
			const double col_j_norm = col_j.norm();
			double temp = fabs(col_i.transpose() * col_j);
			temp *= 1./(col_i_norm*col_j_norm);
			if (mutual_coherence < temp)
				mutual_coherence = temp;
		}
	}
	BOOST_LOG_TRIVIAL(debug)
			<< "Mutual coherence of mapping is " << mutual_coherence;
	return mutual_coherence;
}
