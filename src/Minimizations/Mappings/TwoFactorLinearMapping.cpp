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
 * TwoFactorLinearMapping.cpp
 *
 *  Created on: Mar 16, 2016
 *      Author: heber
 */


#include "BassoConfig.h"

#include "TwoFactorLinearMapping.hpp"

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

TwoFactorLinearMapping::TwoFactorLinearMapping(
		const NormedSpace_weakptr_t _SourceSpaceRef,
		const NormedSpace_weakptr_t _TargetSpaceRef,
		const Eigen::MatrixXd &_first_factor,
		const Eigen::MatrixXd &_second_factor,
		const bool _isAdjoint
		) :
	Mapping(_SourceSpaceRef,_TargetSpaceRef),
	first_factor(_first_factor),
	second_factor(_second_factor),
	isAdjoint(_isAdjoint),
	MatrixVectorProductCounts(0)
{
	assert( second_factor.innerSize() == first_factor.outerSize() );

	if (!isAdjoint) {
		assert( second_factor.outerSize() == getSourceSpace()->getDimension() );
		assert( first_factor.innerSize() == getTargetSpace()->getDimension() );
	} else {
		assert( first_factor.innerSize() == getSourceSpace()->getDimension() );
		assert( second_factor.outerSize() == getTargetSpace()->getDimension() );
	}
}

void TwoFactorLinearMapping::operator()(
		const SpaceElement_ptr_t &_sourceelement,
		SpaceElement_ptr_t &_destelement
		) const
{
	assert( _sourceelement->getSpace() == getSourceSpace() );
	++MatrixVectorProductCounts;

	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();
	Eigen::VectorXd intermediate;
	if (!isAdjoint)
		intermediate = second_factor * RepresentationAdvocate::get(_sourceelement);
	else
		intermediate = first_factor.transpose() * RepresentationAdvocate::get(_sourceelement);
	Eigen::VectorXd tempvector;
	if (!isAdjoint)
		tempvector = first_factor * intermediate;
	else
		tempvector = second_factor.transpose() * intermediate;
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	MatrixVectorProductTimings += timing_end - timing_start;

	RepresentationAdvocate::set(_destelement, tempvector);
}

SpaceElement_ptr_t TwoFactorLinearMapping::operator*(const SpaceElement_ptr_t &_element) const
{
	return operator()(_element);
}

const Mapping_ptr_t TwoFactorLinearMapping::getAdjointMapping() const
{
	if (AdjointTwoFactorLinearMapping.expired()) {
		// create adjoint instance properly
		Mapping_ptr_t adjoint = MappingFactory::createTwoFactorInstance(
				getTargetSpace()->getDualSpace(),
				getSourceSpace()->getDualSpace(),
				first_factor,
				second_factor,
				!isAdjoint);
		const_cast<TwoFactorLinearMapping *>(this)->
				setAdjointMapping(adjoint);
		const_cast<TwoFactorLinearMapping *>(
				static_cast<TwoFactorLinearMapping *>(
						adjoint.get())
						)->setAdjointMapping(SelfRef);
		return adjoint;
	} else {
		// this throws if AdjointTwoFactorLinearMapping is expired
		return Mapping_ptr_t(AdjointTwoFactorLinearMapping);
	}
}

void TwoFactorLinearMapping::setSelfRef(const Mapping_weakptr_t &_selfref)
{
	const_cast<Mapping_weakptr_t &>(SelfRef) = _selfref;
}

void TwoFactorLinearMapping::setAdjointMapping(const Mapping_weakptr_t &_adjoint)
{
	assert( AdjointTwoFactorLinearMapping.expired() );
	const_cast<Mapping_weakptr_t &>(AdjointTwoFactorLinearMapping) = _adjoint;
}

const double TwoFactorLinearMapping::Norm() const
{
#ifdef FULLMATRIXNORM
	Eigen::JacobiSVD<Eigen::MatrixXd> first_svd =
			first_factor.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::JacobiSVD<Eigen::MatrixXd> second_svd =
			second_factor.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);
	const Eigen::JacobiSVD<Eigen::MatrixXd>::SingularValuesType &first_singular_values =
			first_svd.singularValues();
	const Eigen::JacobiSVD<Eigen::MatrixXd>::SingularValuesType &second_singular_values =
			second_svd.singularValues();
	return first_singular_values[0]*second_singular_values[0];
#else
//	if ((matrix.innerSize() != matrix.outerSize())
//			|| (!matrix.isApprox(matrix.transpose()))) {
		LOG(warning, "BEWARE: Is this calculating the right matrix norm?");
//	}
	return first_factor.norm()*second_factor.norm();
#endif
}

const double TwoFactorLinearMapping::MutualCoherence() const
{
	LOG(warning, "BEWARE: We are calculating the matrix product of two factors here for simplicity." << "Do you really need this functionality?");
	const Eigen::MatrixXd matrix = first_factor * second_factor;
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
	LOG(debug, "Mutual coherence of mapping is " << mutual_coherence);
	return mutual_coherence;
}

SingularValueDecomposition TwoFactorLinearMapping::getSVD() const
{
	LOG(warning, "BEWARE: We are calculating the matrix product of two factors here for simplicity." << "Do you really need this functionality?");
	SingularValueDecomposition_impl::ptr_t svd_pimpl(
			new SingularValueDecomposition_impl(first_factor*second_factor));
	SingularValueDecomposition svd(svd_pimpl);
	return svd;
}
