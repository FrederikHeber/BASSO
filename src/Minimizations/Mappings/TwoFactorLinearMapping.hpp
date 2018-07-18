/*
 * TwoFactorLinearMapping.hpp
 *
 *  Created on: Mar 16, 2016
 *      Author: heber
 */

#ifndef TWOFACTORLINEARMAPPING_HPP_
#define TWOFACTORLINEARMAPPING_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>
#include <boost/weak_ptr.hpp>
#include <Eigen/Dense>

#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Mappings/SingularValueDecomposition.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

struct MappingFactory;
class SpaceElement;

/** This class defines a product of two matrices that consecutively
 * transforms elements from one space into another linearly.
 *
 * In detail, given two matrices, the first factor \f$ K \f$ and the second
 * factor \f$ X \f$, we look at \f$ K X x = y \f$ to be considered as
 * \f$ A x = y \f$ in the context of inverse problems, i.e. \f$ A = K X \f$.
 */
class TwoFactorLinearMapping : public Mapping
{
	//!> allow SpaceElement access to matrix
	friend class SpaceElement;

	//!> allow Factory access to private constructors
	friend struct LinearMappingFactory;
private:
	/** Constructor for TwoFactorLinearMapping.
	 *
	 * @param _SourceSpaceRef source space
	 * @param _TargetSpaceRef target space
	 * @param _first_factor first factor \f$ K \f$ of this linear mapping
	 * @param _second_factor second factor \f$ X \f$ of this linear mapping
	 * @param _isAdjoint states whether matrix is to be applied as transposed or not
	 */
	TwoFactorLinearMapping(
			const NormedSpace_weakptr_t _SourceSpaceRef,
			const NormedSpace_weakptr_t _TargetSpaceRef,
			const Eigen::MatrixXd &_first_factor,
			const Eigen::MatrixXd &_second_factor,
			const bool _isAdjoint = false
			);

public:

	virtual ~TwoFactorLinearMapping() {}

	/** Matrix multiplication from the right.
	 *
	 * @param _element element to map/transform
	 * @return new transformed/mapped element
	 */
	SpaceElement_ptr_t operator*(const SpaceElement_ptr_t &_element) const;

	//!> expose overloaded operator method from base class
	using Mapping::operator();

	/** Matrix multiplication from the right.
	 *
	 * @param _sourceelement element to map/transform
	 * @param _destelement transformed/mapped element on return
	 */
	void operator()(
				const SpaceElement_ptr_t &_sourceelement,
				SpaceElement_ptr_t &_destelement
				) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * \note With this specific Mapping (i.e. TwoFactorLinearMapping) we make
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
	{ return MatrixVectorProductCounts; }

	/** Returns the total runtime the program spent so far on
	 * its operator().
	 *
	 * \warning Due to Eigen's lazy evaluation we cannot really
	 * measure the time required on matrix-vector products.
	 *
	 * @return runtime summed over all calls
	 */
	const boost::chrono::nanoseconds getTiming() const
	{ return boost::chrono::nanoseconds(0); }

	/** Calculate the SVD for this linear mapping.
	 *
	 * @return SVD
	 */
	SingularValueDecomposition getSVD() const;

private:
	//!> matrix representation of first factor of this linear mapping
	Eigen::MatrixXd first_factor;

	//!> matrix representation of second factor of this linear mapping
	Eigen::MatrixXd second_factor;

	/** Similar as to internal ref for NormedSpace, the SpaceElements
	 * also holds a weak_ptr reference to itself in order to return
	 * it on specific operators
	 */
	const Mapping_weakptr_t SelfRef;

	//!> internal pointer to adjoint
	const Mapping_weakptr_t AdjointTwoFactorLinearMapping;

	/** Setter for the internal weak_ptr to \a this.
	 *
	 * \note This must be used by the NormedSpace::createElement() as
	 * otherwise certain operators will return shared_ptr's containing NULL.
	 *
	 * @param _selfref reference to this entity wrapped in weak_ptr.
	 */
	void setSelfRef(const Mapping_weakptr_t &_selfref);

	/** Internal setter for the adjoint mapping
	 *
	 * When we create the adjoint mapping for the first time, the
	 * adjoint's adjoint mapping is the original one and hence must
	 * not be created.
	 *
	 * @param _adjoint adjoint instance to set
	 */
	void setAdjointMapping(const Mapping_weakptr_t &_adjoint);

	//!> states whethet this is the normal or the adjoint operator
	const bool isAdjoint;

	//!> typedef for counting mvps
	typedef unsigned int Count_t;

	//!> Counts the number of calculated matrix-vector products
	mutable Count_t MatrixVectorProductCounts;

	//!> counts the time used for matrix-vector products
	mutable boost::chrono::nanoseconds MatrixVectorProductTimings;
};


#endif /* TWOFACTORLINEARMAPPING_HPP_ */
