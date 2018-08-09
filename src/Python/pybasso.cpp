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

#include <boost/python.hpp>
#include <boost/python/operators.hpp>

#include <string>

#include "Minimizations/Elements/SpaceElement.hpp"
#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/DualityMapping.hpp"
#include "Minimizations/Mappings/LinearMapping.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Mappings/NonLinearMapping.hpp"
#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"
#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/types.hpp"
#include "Options/CommandLineOptions.hpp"

#include "Python/utility.hpp"

using namespace boost::python;

/* specific functions for interface */

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

void setVerbosity(CommandLineOptions &opts, unsigned int verbose)
{
	opts.verbose = verbose;
	opts.setVerbosity();
}

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

    /*** classes ***/
    class_<LpNorm, bases<Norm> >("LpNorm", no_init);
    class_<NormedSpace, NormedSpace_ptr_t >("NormedSpace", no_init)
    /* does not work: .def_readonly("norm", &NormedSpace_getNorm) */
        .def("getNorm", &NormedSpace_getNorm, return_internal_reference<1>())
        .def_readonly("space", &NormedSpace::getSpace)
        .def_readonly("dualspace", &NormedSpace::getDualSpace)
        .def("getDualityMapping", &NormedSpace::getDualityMapping)
        .def_readonly("dim", &NormedSpace::getDimension)
        .def("createElement", &NormedSpace::createElement)
    ;

    const bool (SpaceElement::*SpaceElement_isZero_tolerance)(const double) const =
            &SpaceElement::isZero;
    SpaceElement_ptr_t (*double_times_SpaceElement)(
            const double _alpha,
            const SpaceElement_ptr_t &_element) = &operator*;
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
        .def_readonly("space", &SpaceElement::getSpace)
        .def("setZero", &SpaceElement::setZero)
        .def("__getitem__", &SpaceElement_getitem)
        .def("__setitem__", &SpaceElement_setitem)
        .def("get", &pyBasso_SpaceElement_access::get)
        .def("set", &pyBasso_SpaceElement_access::set)
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
    class_<NonLinearMapping, bases<Mapping> >("NonLinearMapping", no_init)
    ;
    class_<CommandLineOptions>("CommandLineOptions", init<>())
        .def_readwrite("algorithm_name", &CommandLineOptions::algorithm_name)
        .def_readwrite("auxiliary_constraints", &CommandLineOptions::auxiliary_constraints)
        .def_readwrite("C", &CommandLineOptions::C)
        .def_readwrite("calculateAngles", &CommandLineOptions::calculateAngles)
        .def_readwrite("delta", &CommandLineOptions::delta)
        .def_readwrite("enforceRandomMapping", &CommandLineOptions::enforceRandomMapping)
        .def_readwrite("everynthtuple", &CommandLineOptions::everynthtuple)
        .def_readwrite("inexactLinesearch", &CommandLineOptions::inexactLinesearch)
        .def_readwrite("iteration_file", &CommandLineOptions::iteration_file)
        .def_readwrite("maxinneriter", &CommandLineOptions::maxinneriter)
        .def_readwrite("maxiter", &CommandLineOptions::maxiter)
        .def_readwrite("max_sfp_loops", &CommandLineOptions::max_sfp_loops)
        .def_readwrite("maxwalltime", &CommandLineOptions::maxwalltime)
        .def_readwrite("minlib", &CommandLineOptions::minlib)
        .def_readwrite("N", &CommandLineOptions::N)
        .add_property("orthogonalization_type",
                      &CommandLineOptions_orthogonalization_type_get, // getter
                      &CommandLineOptions_orthogonalization_type_set) // setter
        .def_readwrite("outputsteps", &CommandLineOptions::outputsteps)
        .def_readwrite("regularization_parameter", &CommandLineOptions::regularization_parameter)
        .def_readwrite("searchspace_type", &CommandLineOptions::searchspace_type)
        .def_readwrite("stepwidth_type", &CommandLineOptions::stepwidth_type)
        .def_readwrite("stopping_criteria", &CommandLineOptions::stopping_criteria)
        .def_readwrite("tau", &CommandLineOptions::tau)
        .def_readwrite("tolerance_linesearch", &CommandLineOptions::tolerance_linesearch)
        .def_readwrite("tolerance_spacex", &CommandLineOptions::tolerance_spacex)
        .def_readwrite("tuple_parameters", &CommandLineOptions::tuple_parameters)
        .def_readwrite("updatetype", &CommandLineOptions::updatetype)
        .def_readwrite("wolfe_constants", &CommandLineOptions::wolfe_constants)
        .def("__repr__", &CommandLineOptions_toString)
        .def("setValues", &CommandLineOptions::setSecondaryValues)
        .def("checkSensibility", &CommandLineOptions::checkSensibility)
        .def("setVerbosity", &setVerbosity)
        ;

    class_<InverseProblem, InverseProblem_ptr_t>("InverseProblem", init<
                const Mapping_ptr_t, const NormedSpace_ptr_t,
                const NormedSpace_ptr_t, const SpaceElement_ptr_t>())
        .def_readonly("sourcespace", &InverseProblem::SourceSpace)
        .def_readonly("targetspace", &InverseProblem::TargetSpace)
        .def_readonly("dualsourcespace", &InverseProblem::DualSourceSpace)
        .def_readonly("dualtargetspace", &InverseProblem::DualTargetSpace)
        .def_readonly("A", &InverseProblem::A)
        .def_readonly("At", &InverseProblem::A_t)
        .def_readwrite("x", &InverseProblem::x)
        .def_readonly("y", &InverseProblem::y)
        .def("solve", &InverseProblem_solve)
        .def("project", &RangeProjection_solve)
    ;

    // factory methods
    def("create_LpSpace", &create_LpSpace);
    def("create_LinearMapping", &create_LinearMapping);
    def("create_NonLinearMapping", &create_NonLinearMapping);
}
