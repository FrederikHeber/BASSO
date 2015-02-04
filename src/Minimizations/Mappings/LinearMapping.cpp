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
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

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
	assert( matrix.outerSize() == getSourceSpace()->getDimension() );
	assert( matrix.innerSize() == getTargetSpace()->getDimension() );
}

const SpaceElement_ptr_t LinearMapping::operator()(
		const SpaceElement_ptr_t &_sourceelement
		) const
{
	assert( _sourceelement->getSpace() == getSourceSpace() );
	SpaceElement_ptr_t targetelement =
			ElementCreator::create(
					getTargetSpace(),
					MatrixVectorProduct(
						matrix,
						 RepresentationAdvocate::get(_sourceelement)));
	return targetelement;
}

SpaceElement_ptr_t LinearMapping::operator*(const SpaceElement_ptr_t &_element) const
{
	assert( _element->getSpace() == getSourceSpace() );
	SpaceElement_ptr_t targetelement =
			ElementCreator::create(
					getTargetSpace(),
					MatrixVectorProduct(
						matrix,
						RepresentationAdvocate::get(_element)));
	return targetelement;
}

const Mapping_ptr_t LinearMapping::getAdjointMapping() const
{
	if (AdjointLinearMapping == NULL) {
		// create adjoint instance properly
		Mapping_ptr_t adjoint = LinearMappingFactory::createInstance(
				getTargetSpace()->getDualSpace(),
				getSourceSpace()->getDualSpace(),
				matrix.transpose());
		const_cast<LinearMapping *>(this)->
				setAdjointMapping(adjoint);
		const_cast<LinearMapping *>(
				static_cast<const LinearMapping *>(
						adjoint.get())
						)->setAdjointMapping(
								Mapping_ptr_t(SelfRef));
	}
	return AdjointLinearMapping;
}

void LinearMapping::setSelfRef(const Mapping_ptr_t &_selfref)
{
	const_cast<boost::weak_ptr<Mapping> &>(SelfRef) = _selfref;
}

void LinearMapping::setAdjointMapping(const Mapping_ptr_t &_adjoint)
{
	assert( AdjointLinearMapping == NULL );
	const_cast<Mapping_ptr_t &>(AdjointLinearMapping) = _adjoint;
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
