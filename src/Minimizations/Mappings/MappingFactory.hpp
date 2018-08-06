/*
 * LinearMappingFactory.hpp
 *
 *  Created on: Dec 12, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MAPPINGS_MAPPINGFACTORY_HPP_
#define MINIMIZATIONS_MAPPINGS_MAPPINGFACTORY_HPP_

#include "BassoConfig.h"

#include <string>

#include <Eigen/Dense>

#include "Minimizations/types.hpp"

/** This factory makes sure that LinearMappings are correctly
 * instantiated.
 *
 * This is required due to the internal LinearMapping::SelfRef
 * which must be set on creation.
 */
struct MappingFactory
{
	/** Factory function creating a linear mapping, i.e. a matrix
	 * associated to a specific space.
	 *
	 * @param _SourceSpaceRef reference to source space
	 * @param _TargetSpaceRef reference to target space
	 * @param _matrix matrix of this linear mapping
	 * @param _isAdjoint states whether matrix is to be applied as transposed or not
	 * @return
	 */
	static const Mapping_ptr_t createInstance(
			const NormedSpace_weakptr_t _SourceSpaceRef,
			const NormedSpace_weakptr_t _TargetSpaceRef,
			const Eigen::MatrixXd & _matrix,
			const bool _isAdjoint = false);

	/** Factory function creating a linear mapping consisting of two matrix
	 * factors whose product is the the linear mapping.
	 *
	 * @param _SourceSpaceRef reference to source space
	 * @param _TargetSpaceRef reference to target space
	 * @param _first_factor first factor of this linear mapping
	 * @param _second_factor second factor of this linear mapping
	 * @param _isAdjoint states whether matrix is to be applied as transposed or not
	 * @return
	 */
	static const Mapping_ptr_t createTwoFactorInstance(
			const NormedSpace_weakptr_t _SourceSpaceRef,
			const NormedSpace_weakptr_t _TargetSpaceRef,
			const Eigen::MatrixXd & _first_factor,
			const Eigen::MatrixXd & _second_factor,
			const bool _isAdjoint = false);

	/** Factory function creating a linear mapping, i.e. a matrix
	 * associated to a specific space.
	 *
	 * @param _SourceSpaceRef reference to source space
	 * @param _TargetSpaceRef reference to target space
	 * @param _name token (or filename) of the matrix object
	 * @param _isAdjoint states whether matrix is to be applied as transposed or not
	 * @return
	 */
	static const Mapping_ptr_t create(
			const NormedSpace_weakptr_t _SourceSpaceRef,
			const NormedSpace_weakptr_t _TargetSpaceRef,
			const std::string & _name,
			const bool _isAdjoint = false);
};



#endif /* MINIMIZATIONS_MAPPINGS_MAPPINGFACTORY_HPP_ */
