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
 * \file utility.cpp
 *
 * In order to deal with boost::weak_ptr issues we have moved all helper
 * functions that use weak_ptrs into this extra module where no boost::python
 * is included. This breaks ODR but helps us to get this part working as it
 * will not use the overridden weak_ptr definitions of boost::python which
 * are not doing the right thing (with respect to the extra reference
 * counting).
 *
 *
 *  Created on: Jul 18, 2018
 *      Author: heber
 */

#include <boost/assign.hpp>

#include "Log/Logging.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblemFactory.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"
#include "Minimizations/types.hpp"

using namespace boost::assign;

NormedSpace_ptr_t create_LpSpace(
		const int _dimension,
		const double _p,
		const double _power
		)
{
	InverseProblemFactory::args_t _args_SpaceX;
	_args_SpaceX += boost::any(_p), boost::any(_power);
	return NormedSpaceFactory::create(
			_dimension, std::string("lp"), _args_SpaceX);
}

const Mapping_ptr_t create_LinearMapping(
		NormedSpace_ptr_t &_SourceSpaceRef,
		NormedSpace_ptr_t &_TargetSpaceRef,
		const Eigen::MatrixXd &_matrix,
		const bool _isAdjoint)
{
	return LinearMappingFactory::createInstance(
			_SourceSpaceRef, _TargetSpaceRef, _matrix, _isAdjoint);
}

double SpaceElement_getitem(const SpaceElement_ptr_t &_element, const int _index)
{ 	return (*_element)[_index]; }

void SpaceElement_setitem(SpaceElement_ptr_t &_element, const int _index, const double _value)
{ 	(*_element)[_index] = _value; }

/* thin wrappers to preserve default arguments */
const bool SpaceElement_isZero(const SpaceElement_ptr_t &_element)
{ return _element->isZero(); }

const bool SpaceElement_isApprox(
		const SpaceElement_ptr_t &_element,
		const SpaceElement_ptr_t &_other,
		const double _tolerance)
{ return _element->isApprox(_other, _tolerance); }

const Norm& NormedSpace_getNorm(const NormedSpace_ptr_t &_space)
{ 	return const_cast<const Norm &>(*_space->getNorm().get()); }

const SpaceElement_ptr_t Mapping_operator(
		const Mapping_ptr_t &_map,
		const SpaceElement_ptr_t &_element)
{ 	return (*_map)(_element); }

const double Mapping_getTiming(const Mapping_ptr_t &_map)
{ 	return boost::chrono::duration<double>(_map->getTiming()).count(); }
