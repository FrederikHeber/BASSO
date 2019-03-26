/*
 * Project: BASSO - BAnach Sequential Subspace Optimizer
 * Description: C++ library for solving optimization problems in Banach spaces
 * Copyright (C)  2014-2016 University of Saarland. All rights reserved.
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
 * MappingFactory.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: heber
 */


#include "MappingFactory.hpp"

#include "BassoConfig.h"

#include <fstream>

#include "MatrixIO/MatrixIO.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/NonLinearMapping.hpp"
#include "Minimizations/Mappings/TwoFactorLinearMapping.hpp"

const Mapping_ptr_t
MappingFactory::createInstance(
		const NormedSpace_weakptr_t _SourceSpaceRef,
		const NormedSpace_weakptr_t _TargetSpaceRef,
		const Eigen::MatrixXd & _matrix,
		const bool _isAdjoint)
{
	Mapping_ptr_t mapping(
			new LinearMapping(
					_SourceSpaceRef,
					_TargetSpaceRef,
					_matrix,
					_isAdjoint));
	static_cast<LinearMapping *>(mapping.get())->setSelfRef(mapping);
	return mapping;
}

const Mapping_ptr_t
MappingFactory::createNonlinearInstance(
		const NormedSpace_weakptr_t _SourceSpaceRef,
		const NormedSpace_weakptr_t _TargetSpaceRef,
		const NonLinearMapping::non_linear_map_t &_map_function,
		const NonLinearMapping::jacobian_t &_derivative,
		const bool _isAdjoint)
{
	Mapping_ptr_t mapping(
			new NonLinearMapping(
					_SourceSpaceRef,
					_TargetSpaceRef,
					_map_function,
					_derivative,
					_isAdjoint));
	static_cast<NonLinearMapping *>(mapping.get())->setSelfRef(mapping);
	return mapping;
}

const Mapping_ptr_t
MappingFactory::createTwoFactorInstance(
		const NormedSpace_weakptr_t _SourceSpaceRef,
		const NormedSpace_weakptr_t _TargetSpaceRef,
		const Eigen::MatrixXd & _first_factor,
		const Eigen::MatrixXd & _second_factor,
		const bool _isAdjoint)
{
	Mapping_ptr_t mapping(
			new TwoFactorLinearMapping(
					_SourceSpaceRef,
					_TargetSpaceRef,
					_first_factor,
					_second_factor,
					_isAdjoint));
	static_cast<TwoFactorLinearMapping *>(mapping.get())->setSelfRef(mapping);
	return mapping;
}

const Mapping_ptr_t
MappingFactory::create(
		const NormedSpace_weakptr_t _SourceSpaceRef,
		const NormedSpace_weakptr_t _TargetSpaceRef,
		const std::string & _name,
		const bool _isAdjoint)
{
	Eigen::MatrixXd matrix;

	using namespace MatrixIO;
	if (MatrixIO::isPresentFile(_name)) {
		std::ifstream ist(_name.c_str());
		if (ist.good()) {
			try {
				ist >> matrix;
			} catch (MatrixIOStreamEnded_exception &e) {
				std::cerr << "Failed to fully parse matrix from " << _name << std::endl;
			}
		} else {
			std::cerr << "Failed to open " << _name << std::endl;
		}
	} else {
		std::cerr << _name << " is not a valid name for SpaceElementFactory." << std::endl;
	}

	return createInstance(
			_SourceSpaceRef,
			_TargetSpaceRef,
			matrix,
			_isAdjoint);
}
