/*
 * LinearMapping.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef LINEARMAPPING_HPP_
#define LINEARMAPPING_HPP_

#include "BassoConfig.h"

#include <cassert>
#include <Eigen/Dense>

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
	/** Default constructor for LinearMapping.
	 *
	 */
	LinearMapping()
	{}

	/** Constructor for LinearMapping.
	 *
	 * @param _SourceSpaceRef source space
	 * @param _TargetSpaceRef target space
	 */
	LinearMapping(
			const NormedSpace_ptr_t _SourceSpaceRef,
			const NormedSpace_ptr_t _TargetSpaceRef
			) :
		Mapping(_SourceSpaceRef,_TargetSpaceRef),
		matrix(Eigen::MatrixXd::Zero(
				_SourceSpaceRef->getDimension(),
				_TargetSpaceRef->getDimension())
		)
	{}

	/** Constructor for LinearMapping with a given matrix
	 *
	 * @param _matrix finite-dimensional representation of the mapping
	 */
	LinearMapping(const Eigen::MatrixXd &_matrix) :
		matrix(_matrix)
	{}

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
			) :
		Mapping(_SourceSpaceRef,_TargetSpaceRef),
		matrix(_matrix)
	{
//		assert( matrix.innerSize() == SourceSpaceRef->getDimension() );
//		assert( matrix.outerSize() == TargetSpaceRef->getDimension() );
	}

	/** Matrix multiplication from the right.
	 *
	 * @param _element element to map/transform
	 * @return new transformed/mapped element
	 */
	SpaceElement_ptr_t operator*(const SpaceElement_ptr_t &_element) const;

	/** Matrix multiplication from the right.
	 *
	 * @param _vector element to map/transform
	 * @return new transformed/mapped vector
	 */
	const Eigen::VectorXd operator*(const Eigen::VectorXd &_vector) const;

	/** Matrix multiplication from the right.
	 *
	 * @param _sourceelement element to map/transform
	 * @return mapped/transformed element
	 */
	SpaceElement_ptr_t operator()(
				const SpaceElement_ptr_t &_sourceelement
				) const;

	/** Mapping function.
	 *
	 * @param _sourcevector vector to map/transform
	 * @return new transformed/mapped vector
	 */
	const Eigen::VectorXd operator()(
			const Eigen::VectorXd &_sourcevector
			) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	Mapping_ptr_t getAdjointMapping() const;

	/** Returns the matrix norm.
	 *
	 * @return norm of \a matrix
	 */
	const double Norm() const;

	/** Const getter for the internal matrix representation.
	 *
	 * @return const ref to matrix
	 */
	const Eigen::MatrixXd& getMatrixRepresentation() const
	{ return matrix; }

private:
	//!> matrix representation of this linear mapping
	Eigen::MatrixXd matrix;
};


#endif /* LINEARMAPPING_HPP_ */
