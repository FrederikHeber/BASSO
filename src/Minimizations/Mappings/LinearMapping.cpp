/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
 * Copyright (C)  2016-2019 Frederik Heber. All rights reserved.
 *
 *
 *   This file is part of the BASSO library.
 *
 *    BASSO is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    BASSO is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with BASSO.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

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
#include "MappingFactory.hpp"

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Mappings/SingularValueDecomposition.hpp"
#include "Minimizations/Mappings/SingularValueDecomposition_impl.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

// static refs
unsigned int LinearMapping::warned_rightnorm = 0;

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
		Mapping_ptr_t adjoint = MappingFactory::createInstance(
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
//			|| (!matrix.isApprox(matrix.transpose()))) {
	if (warned_rightnorm < 3) {
		LOG(warning, "BEWARE: Is this calculating the right matrix norm?");
		++warned_rightnorm;
	}
//	}
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
