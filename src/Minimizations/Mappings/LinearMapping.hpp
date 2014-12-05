/*
 * LinearMapping.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef LINEARMAPPING_HPP_
#define LINEARMAPPING_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "MatrixIO/OperationCounter.hpp"

#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

class SpaceElement;

/** This class defines matrices transforming elements from one
 * space into another linearly.
 */
class LinearMapping : public Mapping
{
	//!> allow SpaceElement access to matrix
	friend class SpaceElement;
public:
	/** Constructor for LinearMapping.
	 *
	 * @param _SourceSpaceRef source space
	 * @param _TargetSpaceRef target space
	 */
	LinearMapping(
			const NormedSpace_ptr_t _SourceSpaceRef,
			const NormedSpace_ptr_t _TargetSpaceRef
			);

	/** Constructor for LinearMapping with a given matrix
	 *
	 * @param _SourceSpaceRef source space
	 * @param _TargetSpaceRef target space
	 * @param _matrix finite-dimensional representation of the mapping
	 */
	LinearMapping(
			const NormedSpace_ptr_t _SourceSpaceRef,
			const NormedSpace_ptr_t _TargetSpaceRef,
			const Eigen::MatrixXd &_matrix
			);

	/** Matrix multiplication from the right.
	 *
	 * @param _element element to map/transform
	 * @return new transformed/mapped element
	 */
	SpaceElement_ptr_t operator*(const SpaceElement_ptr_t &_element) const;

	/** Matrix multiplication from the right.
	 *
	 * @param _sourceelement element to map/transform
	 * @return mapped/transformed element
	 */
	const SpaceElement_ptr_t operator()(
				const SpaceElement_ptr_t &_sourceelement
				) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	const Mapping_ptr_t getAdjointMapping() const;

	/** Returns the matrix norm.
	 *
	 * @return norm of \a matrix
	 */
	const double Norm() const;

	/** Calculate the mutual coherence of the matrix.
	 *
	 * @return maximum of "angle" between column vectors
	 */
	const double MutualCoherence() const;

private:
	//!> matrix representation of this linear mapping
	Eigen::MatrixXd matrix;

	boost::function<
		const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type  (
				const Eigen::MatrixBase<Eigen::MatrixXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&
				)> matrix_vector_fctor;

	const OperationCounter<
		const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type,
		const Eigen::MatrixBase<Eigen::MatrixXd>&,
		const Eigen::MatrixBase<Eigen::VectorXd>&
		> MatrixVectorProduct;
};


#endif /* LINEARMAPPING_HPP_ */
