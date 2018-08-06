/*
 * NonLinearMapping.hpp
 *
 *  Created on: Jul 16, 2018
 *      Author: heber
 */

#ifndef NONLINEARMAPPING_HPP_
#define NONLINEARMAPPING_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>
#include <boost/function.hpp>

#include <Eigen/Dense>

#include "Minimizations/types.hpp"
#include "Minimizations/Mappings/Mapping.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

struct MappingFactory;

/** NonLinearMapping implements a non-linear mapping in the context of inverse
 * problems.
 *
 * Here, we follow the idea of [Wald, Schuster, 2018] to use the transposed
 * derivative of some non-linear operator \f$ F(x) \f$ as the (linearized)
 * adjoint mapping, i.e. we perform updates as \f$ F'(x)^\ast J_p( F(x) -y ) \f$.
 *
 */
class NonLinearMapping : public Mapping
{
	//!> allow Factory access to private constructors
	friend struct MappingFactory;
public:
	//!> typedef for the signature of the mapping function
	typedef boost::function<Eigen::VectorXd (const Eigen::VectorXd &)> non_linear_map_t;

	/** Constructor for a non-linear mapping.
	 *
	 * @param _SourceSpaceRef source space reference
	 * @param _TargetSpaceRef target space reference
	 */
	NonLinearMapping(
			const NormedSpace_weakptr_t &_SourceSpaceRef,
			const NormedSpace_weakptr_t &_TargetSpaceRef,
			const non_linear_map_t &_map_function,
			const non_linear_map_t &_derivative,
			const bool _isAdjoint
			);

	virtual ~NonLinearMapping() {}

	/** Non-linear mapping function.
	 *
	 * @param _sourceelement element to map/transform
	 * @param _destelement transformed/mapped element on return
	 */
	virtual void operator()(
			const SpaceElement_ptr_t &_sourceelement,
			SpaceElement_ptr_t &_destelement
			) const;

	/** Returns the transposed derivative of the non-linear mapping as the
	 * linearized adjoint of the operator.
	 *
	 * @return transposed derivative of the non-linear mapping
	 */
	virtual const Mapping_ptr_t getAdjointMapping() const;

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

protected:
	//!> reference to space from which this mappings projects
	const NormedSpace_weakptr_t SourceSpaceRef;

	//!> reference to space into which this mappings projects
	const NormedSpace_weakptr_t TargetSpaceRef;

	//!> reference to the non-linear map function
	const non_linear_map_t map_function;

	//!> reference to the derivative of the non-linear map
	const non_linear_map_t derivative;

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



#endif /* NONLINEARMAPPING_HPP_ */
