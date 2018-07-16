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

/* pybasso.cpp
 *
 *  Created on: Jul 16, 2018
 *      Author: heber
 */

#include <boost/assign.hpp>
#include <boost/python.hpp>
#include <boost/python/operators.hpp>

#include <string>

#include "Log/Logging.hpp"

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblemFactory.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/LinearMappingFactory.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/NormedSpaceFactory.hpp"
#include "Minimizations/types.hpp"

using namespace boost::assign;
using namespace boost::python;

/* specific functions for interface */

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
		NormedSpace_const_ptr_t &_SourceSpaceRef,
		NormedSpace_const_ptr_t &_TargetSpaceRef,
		const Eigen::MatrixXd &_matrix,
		const bool _isAdjoint)
{
	LOG(info, "Source space's weak ptr has " << _SourceSpaceRef.use_count());
	LOG(info, "Source space is located at " << NormedSpace_ptr_t(_SourceSpaceRef).get());
	LOG(info, "Target space is located at " << NormedSpace_ptr_t(_TargetSpaceRef).get());
	LOG(info, "Dual source space is located at " << NormedSpace_ptr_t(_SourceSpaceRef)->getDualSpace().get());
	LOG(info, "Dual target space is located at " << NormedSpace_ptr_t(_TargetSpaceRef)->getDualSpace().get());
	Eigen::MatrixXd matrix(10,10);
	return LinearMappingFactory::createInstance(
			_SourceSpaceRef->shared_from_this(), _TargetSpaceRef->shared_from_this(), _matrix, _isAdjoint);
}

double SpaceElement_getitem(const SpaceElement_ptr_t &_element, const int _index)
{
	return (*_element)[_index];
}

void SpaceElement_setitem(SpaceElement_ptr_t &_element, const int _index, const double _value)
{
	(*_element)[_index] = _value;
}

/* thin wrappers to preserve default arguments */
const bool SpaceElement_isZero(const SpaceElement_ptr_t &_element) { return _element->isZero(); }

const bool SpaceElement_isApprox(
		const SpaceElement_ptr_t &_element,
		const SpaceElement_ptr_t &_other,
		const double _tolerance) { return _element->isApprox(_other, _tolerance); }

const Norm& NormedSpace_getNorm(const NormedSpace_ptr_t &_space) {
	return const_cast<const Norm &>(*_space->getNorm().get()); }

/* Wrapper class for interface class with virtual functions */

struct NormWrap : Norm, wrapper<Norm> {
    virtual bool isSmooth() const
    { return this->get_override("isSmooth")(); }
    virtual const double internal_operator(const SpaceElement_ptr_t &_element) const
    { return this->get_override("internal_operator")(); }
};

struct MappingWrap : Mapping, wrapper<Mapping> {
	virtual void operator()(
				const SpaceElement_ptr_t &_sourceelement,
				SpaceElement_ptr_t &_destelement
				)
    { this->get_override("operator()")(); }
	virtual const Mapping_ptr_t getAdjointMapping()
    { return this->get_override("getAdjointMapping")(); }
	virtual const unsigned int getCount()
    { return this->get_override("getCount")(); }
	virtual const boost::chrono::nanoseconds getTiming()
    { return this->get_override("getTiming")(); }
};

BOOST_PYTHON_MODULE(pyBasso)
{
	/*** classes with virtual functions ***/
	class_<NormWrap, Norm_ptr_t, boost::noncopyable>("Norm", no_init)
		.def_readonly("p", &Norm::getPvalue)
		.def_readonly("space", &Norm::getSpace)
		.def("__call__", &Norm::operator())
	    ;
	class_<MappingWrap, Mapping_ptr_t, boost::noncopyable>("Mapping", no_init)
		.def_readonly("sourcespace", &Mapping::getSourceSpace)
		.def_readonly("targetspace", &Mapping::getTargetSpace)
		.def("getSourceSpace", &Mapping::getSourceSpace)
/*		.def("__call__", &Mapping::operator()) */
		.def_readonly("adjoint", pure_virtual(&Mapping::getAdjointMapping)) //, return_internal_reference<1>()
		.def_readonly("power", &Mapping::getPower)
		.def("getCount", pure_virtual(&Mapping::getAdjointMapping))
		.def("getTiming", pure_virtual(&Mapping::getAdjointMapping))
	    ;

	/*** classes ***/
	class_<LpNorm, bases<Norm> >("LpNorm", no_init);
    class_<NormedSpace, NormedSpace_ptr_t >("NormedSpace", no_init)
	/* does not work: .def_readonly("norm", &NormedSpace_getNorm) */
		.def("getNorm", &NormedSpace_getNorm, return_internal_reference<1>())
		.def_readonly("space", &NormedSpace::getSpace)
		.def_readonly("dualspace", &NormedSpace::getDualSpace)
		.def("getDualityMapping", &NormedSpace::getDualityMapping, return_internal_reference<1>())
		.def_readonly("dim", &NormedSpace::getDimension)
		.def("createElement", &NormedSpace::createElement)
    ;

    const bool (SpaceElement::*SpaceElement_isZero_tolerance)(const double) const =
    		&SpaceElement::isZero;
    class_<SpaceElement, SpaceElement_ptr_t >("SpaceElement", no_init)
        .def("isInSpace", &SpaceElement::isInSpace)
        .def("isZero", SpaceElement_isZero)
        .def("isZero", SpaceElement_isZero_tolerance)
        .def("isApproxToConstant", &SpaceElement::isApproxToConstant)
        .def("isApprox", &SpaceElement_isApprox)
        .def("isNonnegative", &SpaceElement::isNonnegative)
        .def("getSignVector", &SpaceElement::getSignVector)
        .def("getAbsVector", &SpaceElement::getAbsVector)
        .def("getCircShiftedVector", &SpaceElement::getCircShiftedVector)
        .def("getMaxCoefficientAndIndex", &SpaceElement::getMaxCoefficientAndIndex)
        .def_readonly("norm", &SpaceElement::Norm)
        .def("pow", &SpaceElement::pow)
        .def("getSpace", &SpaceElement::getSpace)
        .def("setZero", &SpaceElement::setZero)
		.def("__getitem__", &SpaceElement_getitem)
		.def("__setitem__", &SpaceElement_setitem)
		.def(self_ns::str(self_ns::self))
		.def(self_ns::self * float())
		.def(self_ns::self + SpaceElement_ptr_t())
		.def(self_ns::self += SpaceElement_ptr_t())
		.def(self_ns::self - SpaceElement_ptr_t())
		.def(self_ns::self -= SpaceElement_ptr_t())
		.def(self_ns::self == self_ns::self)
	;
	class_<LinearMapping, bases<Mapping> >("LinearMapping", no_init)
		.def(self_ns::self * SpaceElement_ptr_t())
	;

    // factory methods
    def("create_LpSpace", &create_LpSpace);
    def("create_LinearMapping", &create_LinearMapping);
}
