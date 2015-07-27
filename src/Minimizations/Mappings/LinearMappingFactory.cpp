/*
 * LinearMappingFactory.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "LinearMappingFactory.hpp"

#include <fstream>

#include "MatrixIO/MatrixIO.hpp"
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

const Mapping_ptr_t
LinearMappingFactory::create(
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
