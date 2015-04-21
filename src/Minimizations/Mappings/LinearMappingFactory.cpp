/*
 * LinearMappingFactory.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "LinearMappingFactory.hpp"

#include "Minimizations/Mappings/LinearMapping.hpp"

const Mapping_ptr_t
LinearMappingFactory::createInstance(
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

