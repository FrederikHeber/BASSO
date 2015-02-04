/*
 * LinearMapping.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef LINEARMAPPING_HPP_
#define LINEARMAPPING_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>
#include <boost/weak_ptr.hpp>
#include <Eigen/Dense>

#include "MatrixIO/OperationCounter.hpp"

#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

struct LinearMappingFactory;
class SpaceElement;

/** This class defines matrices transforming elements from one
 * space into another linearly.
 */
class LinearMapping : public Mapping
{
	//!> allow SpaceElement access to matrix
	friend class SpaceElement;

	//!> allow Factory access to private constructors
	friend struct LinearMappingFactory;
private:
	/** Constructor for LinearMapping.
	 *
	 * @param _SourceSpaceRef source space
	 * @param _TargetSpaceRef target space
	 * @param _matrix finite-dimensional representation of the mapping
	 */
	LinearMapping(
			const NormedSpace_weakptr_t _SourceSpaceRef,
			const NormedSpace_weakptr_t _TargetSpaceRef,
			const Eigen::MatrixXd &_matrix
			);

public:

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
	 * \note With this specific Mapping (i.e. LinearMapping) we make
	 * sure that only a single adjoint mapping (per individual mapping)
	 * is ever created.
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

	/** Returns the number of times the operator() was called.
	 *
	 * @return number of calls
	 */
	const unsigned int getCount() const
	{ return MatrixVectorProduct.getCount(); }

	/** Returns the total runtime the program spent so far on
	 * its operator().
	 *
	 * @return runtime summed over all calls
	 */
	const boost::chrono::nanoseconds getTiming() const
	{ return MatrixVectorProduct.getTiming(); }

private:
	//!> matrix representation of this linear mapping
	Eigen::MatrixXd matrix;

	/** Similar as to internal ref for NormedSpace, the SpaceElements
	 * also holds a weak_ptr reference to itself in order to return
	 * it on specific operators
	 */
	const Mapping_weakptr_t SelfRef;

	//!> internal pointer to adjoint
	const Mapping_weakptr_t AdjointLinearMapping;

	/** Setter for the internal weak_ptr to \a this.
	 *
	 * \note This must be used by the NormedSpace::createElement() as
	 * otherwise certain operators will return shared_ptr's containing NULL.
	 *
	 * @param _selfref reference to this entity wrapped in weak_ptr.
	 */
	void setSelfRef(const Mapping_ptr_t &_selfref);

	/** Internal setter for the adjoint mapping
	 *
	 * When we create the adjoint mapping for the first time, the
	 * adjoint's adjoint mapping is the original one and hence must
	 * not be created.
	 *
	 * @param _adjoint adjoint instance to set
	 */
	void setAdjointMapping(const Mapping_weakptr_t &_adjoint);

	//!> internally bound function to count matrix vector products
	boost::function<
		const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type  (
				const Eigen::MatrixBase<Eigen::MatrixXd>&,
				const Eigen::MatrixBase<Eigen::VectorXd>&
				)> matrix_vector_fctor;

public:
	//!> matrix vector product counter
	const OperationCounter<
		const Eigen::ProductReturnType<Eigen::MatrixXd, Eigen::VectorXd>::Type,
		const Eigen::MatrixBase<Eigen::MatrixXd>&,
		const Eigen::MatrixBase<Eigen::VectorXd>&
		> MatrixVectorProduct;
};


#endif /* LINEARMAPPING_HPP_ */
