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

LinearMapping::LinearMapping(
		const NormedSpace_ptr_t _SourceSpaceRef,
		const NormedSpace_ptr_t _TargetSpaceRef
		) :
	Mapping(_SourceSpaceRef,_TargetSpaceRef),
	matrix(Eigen::MatrixXd::Zero(
			_SourceSpaceRef->getDimension(),
			_TargetSpaceRef->getDimension())
	),
	matrix_vector_fctor(
			boost::bind(
					static_cast<const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type
						(Eigen::MatrixBase<Eigen::MatrixXd>::*)(const Eigen::MatrixBase<Eigen::VectorXd>&) const>(
								&Eigen::MatrixBase<Eigen::MatrixXd>::operator*),
								_1, _2
			)
	),
	MatrixVectorProduct(matrix_vector_fctor)
{}

LinearMapping::LinearMapping(
		const NormedSpace_ptr_t _SourceSpaceRef,
		const NormedSpace_ptr_t _TargetSpaceRef,
		const Eigen::MatrixXd &_matrix
		) :
	Mapping(_SourceSpaceRef,_TargetSpaceRef),
	matrix(_matrix),
	matrix_vector_fctor(
			boost::bind(
					static_cast<const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type
						(Eigen::MatrixBase<Eigen::MatrixXd>::*)(const Eigen::MatrixBase<Eigen::VectorXd>&) const>(
								&Eigen::MatrixBase<Eigen::MatrixXd>::operator*),
								_1, _2
			)
	),
	MatrixVectorProduct(matrix_vector_fctor)
{
	assert( matrix.outerSize() == SourceSpaceRef->getDimension() );
	assert( matrix.innerSize() == TargetSpaceRef->getDimension() );
}

const SpaceElement_ptr_t LinearMapping::operator()(
		const SpaceElement_ptr_t &_sourceelement
		) const
{
	assert( _sourceelement->getSpace() == SourceSpaceRef );
	SpaceElement_ptr_t targetelement =
			ElementCreator::create(
					TargetSpaceRef,
					MatrixVectorProduct(
						matrix,
						 RepresentationAdvocate::get(_sourceelement)));
	return targetelement;
}

SpaceElement_ptr_t LinearMapping::operator*(const SpaceElement_ptr_t &_element) const
{
	assert( _element->getSpace() == SourceSpaceRef );
	SpaceElement_ptr_t targetelement =
			ElementCreator::create(
					TargetSpaceRef,
					MatrixVectorProduct(
						matrix,
						RepresentationAdvocate::get(_element)));
	return targetelement;
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
