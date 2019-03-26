/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2018 Frederik Heber. All rights reserved.
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
 * NonLinearMapping.cpp
 *
 *  Created on: Jul 16, 2018
 *      Author: heber
 */


#include "BassoConfig.h"

#include "NonLinearMapping.hpp"

#include <cassert>
#include <Eigen/Dense>

#include "Log/Logging.hpp"
#include "MappingFactory.hpp"

#include "Minimizations/Elements/ElementCreator.hpp"
#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/Elements/RepresentationAdvocate.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

NonLinearMapping::NonLinearMapping(
		const NormedSpace_weakptr_t &_SourceSpaceRef,
		const NormedSpace_weakptr_t &_TargetSpaceRef,
		const non_linear_map_t &_map_function,
		const jacobian_t &_derivative,
		const bool _isAdjoint
		) :
	Mapping(_SourceSpaceRef,_TargetSpaceRef),
	map_function(_map_function),
	derivative(_derivative),
	isAdjoint(_isAdjoint),
	MatrixVectorProductCounts(0)
{
	// the adjoint of a non-linear mapping is a (approximate) linear mapping
	assert( !isAdjoint );
}

void NonLinearMapping::operator()(
		const SpaceElement_ptr_t &_sourceelement,
		SpaceElement_ptr_t &_destelement
		) const
{
	assert( _sourceelement->getSpace() == getSourceSpace() );
	++MatrixVectorProductCounts;

	const boost::chrono::high_resolution_clock::time_point timing_start =
			boost::chrono::high_resolution_clock::now();

	const Eigen::VectorXd &sourcevector = RepresentationAdvocate::get(_sourceelement);
	const Eigen::VectorXd tempvector =
			map_function( sourcevector );

	const boost::chrono::high_resolution_clock::time_point timing_end =
			boost::chrono::high_resolution_clock::now();
	MatrixVectorProductTimings += timing_end - timing_start;

	RepresentationAdvocate::set(_destelement, tempvector);
}

void NonLinearMapping::updateAdjoint(
		const SpaceElement_ptr_t &_element
		) const
{
	// calculate gradient
	if (!AdjointLinearMapping.expired()) {
		Mapping_ptr_t adjoint = Mapping_ptr_t(AdjointLinearMapping);
		const Eigen::MatrixXd jacobian = derivative(
				RepresentationAdvocate::get(_element) );
		static_cast<LinearMapping *>(adjoint.get())->matrix = jacobian;
	}
}

const Mapping_ptr_t NonLinearMapping::getAdjointMapping() const
{
	if (AdjointLinearMapping.expired()) {
		// use zero as initial jacobian
		const Eigen::MatrixXd jacobian(
				getSourceSpace()->getDualSpace()->getDimension(),
				getTargetSpace()->getDualSpace()->getDimension()
				);
		// create adjoint instance as LinearMapping with the jacobian
		Mapping_ptr_t derivative_as_adjoint = Mapping_ptr_t(
				new LinearMapping(
						getTargetSpace()->getDualSpace(),
						getSourceSpace()->getDualSpace(),
						jacobian,
						false /* this just sets the transpose */));
		static_cast<LinearMapping *>(derivative_as_adjoint.get())->
				setSelfRef(derivative_as_adjoint);

		const_cast<NonLinearMapping *>(this)->
				setAdjointMapping(derivative_as_adjoint);
		const_cast<LinearMapping *>(
				static_cast<LinearMapping *>(
						derivative_as_adjoint.get())
						)->setAdjointMapping(SelfRef);

		// update adjoint
		updateAdjoint(getSourceSpace()->createElement());

		return derivative_as_adjoint;
	} else {
		// this throws if AdjointLinearMapping is expired
		return Mapping_ptr_t(AdjointLinearMapping);
	}
}

void NonLinearMapping::setSelfRef(const Mapping_weakptr_t &_selfref)
{
	const_cast<Mapping_weakptr_t &>(SelfRef) = _selfref;
}

void NonLinearMapping::setAdjointMapping(const Mapping_weakptr_t &_adjoint)
{
	assert( AdjointLinearMapping.expired() );
	const_cast<Mapping_weakptr_t &>(AdjointLinearMapping) = _adjoint;
}
