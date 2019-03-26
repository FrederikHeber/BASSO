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
 * Backprojection.cpp
 *
 *  Created on: Jul 9, 2015
 *      Author: heber
 */


#include "BassoConfig.h"

#include "Backprojection.hpp"

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"

Backprojection::Backprojection(
		const NormedSpace_weakptr_t _SourceSpaceRef,
		const NormedSpace_weakptr_t _TargetSpaceRef,
		const unsigned int _num_pixel_x,
		const unsigned int _num_pixel_y,
		const unsigned int _num_angles,
		const unsigned int _num_offsets) :
	Mapping(_SourceSpaceRef, _TargetSpaceRef),
	backprojection_matrix(_num_pixel_x, _num_pixel_y, _num_angles, _num_offsets),
	MatrixVectorProductCounts(0)
{}

void Backprojection::operator()(
		const SpaceElement_ptr_t &_sourceelement,
		SpaceElement_ptr_t &_destelement
		) const
{
	assert( _sourceelement->getSpace() == getSourceSpace() );
	++MatrixVectorProductCounts;

	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();
	Eigen::VectorXd tempvector;
	tempvector = backprojection_matrix.getMatrix() * RepresentationAdvocate::get(_sourceelement);
	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	MatrixVectorProductTimings += timing_end - timing_start;

	RepresentationAdvocate::set(_destelement, tempvector);
}

const Mapping_ptr_t Backprojection::getAdjointMapping() const
{
	// this throws if AdjointLinearMapping is expired
	return Mapping_ptr_t(AdjointLinearMapping);
}

void Backprojection::setAdjointMapping(const Mapping_weakptr_t &_adjoint)
{
	assert( AdjointLinearMapping.expired() );
	const_cast<Mapping_weakptr_t &>(AdjointLinearMapping) = _adjoint;
}
