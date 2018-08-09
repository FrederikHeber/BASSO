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
 * Mapping.cpp
 *
 *  Created on: Aug 9, 2018
 *      Author: heber
 */

#include <boost/python.hpp>

#include "Minimizations/Mappings/DualityMapping.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Mappings/NonLinearMapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"

#include "Python/utility.hpp"

using namespace boost::python;

/** We need to tell boost::python how to convert the python function
 * into the correct boost::function wrapped function and we need to
 * force this conversion to occur on the python side (using the
 * minieigen module through extract<>()) and not on the C++ side that
 * does not know how to do it.
 *
 * This is taken from https://stackoverflow.com/a/2181710/1967646
 */
struct non_linear_map_wrapper_t
{
	non_linear_map_wrapper_t( object callable ) : _callable( callable ) {}

	Eigen::VectorXd operator()(const Eigen::VectorXd & _arg)
    {
        // These GIL calls make it thread safe, may or may not be needed depending on your use case
        PyGILState_STATE gstate = PyGILState_Ensure();
        Eigen::VectorXd ret = extract<Eigen::VectorXd>(_callable(_arg));
        PyGILState_Release( gstate );
        return ret;
    }

    object _callable;
};

NonLinearMapping::non_linear_map_t createNonLinearMapWrapper( object function )
{
  return NonLinearMapping::non_linear_map_t( non_linear_map_wrapper_t( function ) );
}

struct jacobian_wrapper_t
{
	jacobian_wrapper_t( object callable ) : _callable( callable ) {}

	Eigen::MatrixXd operator()(const Eigen::VectorXd & _arg)
    {
        // These GIL calls make it thread safe, may or may not be needed depending on your use case
        PyGILState_STATE gstate = PyGILState_Ensure();
        Eigen::MatrixXd ret = extract<Eigen::MatrixXd>(_callable(_arg));
        PyGILState_Release( gstate );
        return ret;
    }

    object _callable;
};

NonLinearMapping::jacobian_t createJacobianWrapper( object function )
{
  return NonLinearMapping::jacobian_t( jacobian_wrapper_t( function ) );
}

// make sure that numpy -> Eigen conversion is happening in python code
const Mapping_ptr_t create_NonLinearMapping(
		NormedSpace_ptr_t &_SourceSpaceRef,
		NormedSpace_ptr_t &_TargetSpaceRef,
		const object &_map_function,
		const object &_derivative,
		const bool _isAdjoint)
{
	return create_NonLinearMapping_full(
			_SourceSpaceRef,
			_TargetSpaceRef,
			createNonLinearMapWrapper(_map_function),
			createJacobianWrapper(_derivative),
			_isAdjoint);
}

/* Wrapper class for interface class with virtual functions */

struct MappingWrap : Mapping, wrapper<Mapping> {
    virtual void operator()(
                const SpaceElement_ptr_t &_sourceelement,
                SpaceElement_ptr_t &_destelement
                )
    { this->get_override("operator()")(); }
    virtual const NormedSpace_ptr_t getSourceSpace() const
    { return this->get_override("getSourceSpace")(); }
    virtual const NormedSpace_ptr_t getTargetSpace() const
    { return this->get_override("getTargetSpace")(); }
    virtual const Mapping_ptr_t getAdjointMapping()
    { return this->get_override("getAdjointMapping")(); }
    virtual const unsigned int getCount()
    { return this->get_override("getCount")(); }
    virtual const boost::chrono::nanoseconds getTiming()
    { return this->get_override("getTiming")(); }
};

struct DualityMappingWrap : DualityMapping, wrapper<DualityMapping> {
    virtual void operator()(
                const SpaceElement_ptr_t &_sourceelement,
                SpaceElement_ptr_t &_destelement
                )
    { this->get_override("operator()")(); }
};

struct PowerTypeDualityMappingWrap : PowerTypeDualityMapping, wrapper<PowerTypeDualityMapping> {
    virtual void operator()(
                const SpaceElement_ptr_t &_sourceelement,
                SpaceElement_ptr_t &_destelement
                )
    { this->get_override("operator()")(); }
};


void export_mapping()
{
	class_<MappingWrap, Mapping_ptr_t, boost::noncopyable>("Mapping", no_init)
        .def_readonly("sourcespace", &Mapping::getSourceSpace)
        .def_readonly("targetspace", &Mapping::getTargetSpace)
        .def("__call__", &Mapping_operator)
        .def_readonly("adjoint", &Mapping::getAdjointMapping)
        .def_readonly("power", &Mapping::getPower)
        .def("getCount", &Mapping::getCount)
        .def("getTiming", &Mapping_getTiming)
        ;
    class_<
        DualityMappingWrap,
        boost::shared_ptr<DualityMappingWrap>,
        boost::noncopyable>("DualityMapping", no_init)
        ;
    class_<
        PowerTypeDualityMappingWrap,
        boost::shared_ptr<PowerTypeDualityMappingWrap>,
        boost::noncopyable>("PowerTypeDualityMapping", no_init)
        ;

    class_<LinearMapping, bases<Mapping> >("LinearMapping", no_init)
        .def(self_ns::self * SpaceElement_ptr_t())
    ;
    class_<NonLinearMapping, bases<Mapping> >("NonLinearMapping", no_init)
    ;

    def("create_LinearMapping", &create_LinearMapping);
    def("create_NonLinearMapping", &create_NonLinearMapping);
}
