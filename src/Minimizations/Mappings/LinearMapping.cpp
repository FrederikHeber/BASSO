/*
 * LinearMapping.cpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "LinearMapping.hpp"

#include <cassert>
#include <Eigen/Dense>

#include "Log/Logging.hpp"

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/Mappings/SingularValueDecomposition.hpp"
#include "Minimizations/Mappings/SingularValueDecomposition_impl.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

LinearMapping::LinearMapping(
		const NormedSpace_weakptr_t _SourceSpaceRef,
		const NormedSpace_weakptr_t _TargetSpaceRef,
		const Eigen::MatrixXd &_matrix,
		const bool _isAdjoint
		) :
	Mapping(_SourceSpaceRef,_TargetSpaceRef),
	matrix(_matrix),
	isAdjoint(_isAdjoint),
	MatrixVectorProductCounts(0)
{
	if (!isAdjoint) {
		assert( matrix.outerSize() == getSourceSpace()->getDimension() );
		assert( matrix.innerSize() == getTargetSpace()->getDimension() );
	} else {
		assert( matrix.innerSize() == getSourceSpace()->getDimension() );
		assert( matrix.outerSize() == getTargetSpace()->getDimension() );
	}
}

void LinearMapping::operator()(
		const SpaceElement_ptr_t &_sourceelement,
		SpaceElement_ptr_t &_destelement
		) const
{
	assert( _sourceelement->getSpace() == getSourceSpace() );
	++MatrixVectorProductCounts;

	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();
	Eigen::VectorXd tempvector;
	if (!isAdjoint)
		tempvector = matrix * RepresentationAdvocate::get(_sourceelement);
	else
		tempvector = matrix.transpose() * RepresentationAdvocate::get(_sourceelement);
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	MatrixVectorProductTimings += timing_end - timing_start;

	RepresentationAdvocate::set(_destelement, tempvector);
}

SpaceElement_ptr_t LinearMapping::operator*(const SpaceElement_ptr_t &_element) const
{
	return operator()(_element);
}

const Mapping_ptr_t LinearMapping::getAdjointMapping() const
{
	if (AdjointLinearMapping.expired()) {
		// create adjoint instance properly
		Mapping_ptr_t adjoint = LinearMappingFactory::createInstance(
				getTargetSpace()->getDualSpace(),
				getSourceSpace()->getDualSpace(),
				matrix,
				!isAdjoint);
		const_cast<LinearMapping *>(this)->
				setAdjointMapping(adjoint);
		const_cast<LinearMapping *>(
				static_cast<LinearMapping *>(
						adjoint.get())
						)->setAdjointMapping(SelfRef);
		return adjoint;
	} else {
		// this throws if AdjointLinearMapping is expired
		return Mapping_ptr_t(AdjointLinearMapping);
	}
}

void LinearMapping::setSelfRef(const Mapping_weakptr_t &_selfref)
{
	const_cast<Mapping_weakptr_t &>(SelfRef) = _selfref;
}

void LinearMapping::setAdjointMapping(const Mapping_weakptr_t &_adjoint)
{
	assert( AdjointLinearMapping.expired() );
	const_cast<Mapping_weakptr_t &>(AdjointLinearMapping) = _adjoint;
}

const double LinearMapping::Norm() const
{
#ifdef FULLMATRIXNORM
	Eigen::JacobiSVD<Eigen::MatrixXd> svd =
			matrix.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
	const Eigen::JacobiSVD<Eigen::MatrixXd>::SingularValuesType &singular_values =
			svd.singularValues();
	return singular_values[0];
#else
//	if ((matrix.innerSize() != matrix.outerSize())
//			|| (!matrix.isApprox(matrix.transpose())))
		LOG(warning, "BEWARE: Is this calculating the right matrix norm?");
	return matrix.norm();
#endif
}

const double LinearMapping::MutualCoherence() const
{
	double mutual_coherence = 0.;
	const unsigned int dim = getSourceSpace()->getDimension();
	for (unsigned int i=0;i<dim;++i) {
		for (unsigned int j=i+1;j<dim;++j) {
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
	LOG(debug, "Mutual coherence of mapping is " << mutual_coherence);
	return mutual_coherence;
}

SingularValueDecomposition LinearMapping::getSVD() const
{
	SingularValueDecomposition_impl::ptr_t svd_pimpl(
			new SingularValueDecomposition_impl(matrix));
	SingularValueDecomposition svd(svd_pimpl);
	return svd;
}
