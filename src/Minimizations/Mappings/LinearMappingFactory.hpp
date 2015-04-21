/*
 * LinearMappingFactory.hpp
 *
 *  Created on: Dec 12, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_MAPPINGS_LINEARMAPPINGFACTORY_HPP_
#define MINIMIZATIONS_MAPPINGS_LINEARMAPPINGFACTORY_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/types.hpp"

/** This factory makes sure that LinearMappings are correctly
 * instantiated.
 *
 * This is required due to the internal LinearMapping::SelfRef
 * which must be set on creation.
 */
struct LinearMappingFactory
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
};



#endif /* MINIMIZATIONS_MAPPINGS_LINEARMAPPINGFACTORY_HPP_ */
