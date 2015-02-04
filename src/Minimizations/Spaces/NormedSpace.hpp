/*
 * NormedSpace.hpp
 *
 *  Created on: Oct 23, 2014
 *      Author: heber
 */

#ifndef NORMEDSPACE_HPP_
#define NORMEDSPACE_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

#include <boost/function.hpp>
#include <boost/weak_ptr.hpp>
#include <Eigen/Dense>

#include "MatrixIO/OperationCounter.hpp"

#include "Minimizations/Spaces/VectorSpaceOperationCounts.hpp"

class Norm;
class NormedSpaceFactory;
class SpaceElement;

class NormedSpace
{
private:
	//!> grant factory access to default cstor.
	friend class NormedSpaceFactory;

	//!> typedef for the (lazy) dualitymapping constructor
	typedef boost::function<Mapping_ptr_t()> constructDualityMapping_t;

	/** Private default constructor to which only NormedSpaceFactory has access.
	 *
	 * @param _dimension dimension of this space
	 */
	NormedSpace(
			const unsigned int _dimension);

	/** Setter of the internal space ref.
	 *
	 * \note This is necessary due to current interface limitation working
	 * with shared_ptr's. We must not create a shared_ptr with \a this as
	 * this shared_ptr already exists and must be supplied elsewhere.
	 *
	 * \warning The NormedSpaceFactory function should make sure that a correct
	 * space (with the shared_ptr ref) has been created.
	 *
	 * @param _space space reference
	 */
	void setSpace(const NormedSpace_ptr_t& _space)
	{	const_cast<boost::weak_ptr<NormedSpace> &>(Space) = _space; }

	/** Setter of the dual space.
	 *
	 * This is for the factory only to allow it to connect the various
	 * spaces among one another.
	 *
	 * @param _dualspace
	 */
	void setDualSpace(const NormedSpace_ptr_t& _dualspace)
	{	const_cast<NormedSpace_ptr_t&>(DualSpace) = _dualspace;	}

	/** Setter of the norm.
	 *
	 * This is for the factory only to allow it to connect the various
	 * spaces among one another.
	 *
	 * @param _norm norm object
	 */
	void setNorm(const Norm_ptr_t& _norm)
	{	const_cast<Norm_ptr_t &>(norm) = _norm;	}

	/** Setter of the duality mapping.
	 *
	 * This is for the factory only to allow it to connect the various
	 * spaces among one another.
	 *
	 * @param _mapping duality mapping object
	 */
	void setDualityMapping(const Mapping_ptr_t &_mapping)
	{	const_cast<Mapping_ptr_t &>(dualitymapping) = _mapping;	}

	/** Setter of the duality mapping constructor.
	 *
	 * We need the constructor to overcome the problem of space and dual
	 * space not being constructed at the same time.
	 *
	 * @param _constructDualityMapping duality mapping object constructor
	 */
	void setDualityMappingConstructor(
			const constructDualityMapping_t &_constructDualityMapping)
	{	const_cast<constructDualityMapping_t &>(constructDualityMapping) =
			_constructDualityMapping;
	}

public:
	/** Constructor of NormedSpace.
	 *
	 * \warning setSpace() must be called to provide shared_ptr instace of
	 * this object. Otherwise getSpace() will return null and 
	 * createElement() does not work, rendering this space unuseable.
	 *
	 * @param _dimension dimension of this space
	 * @param _norm norm object of the space
	 * @param _constructDualityMapping stored function call for constructing
	 * 		\a dualitymapping on demand
	 */
	NormedSpace(
			const unsigned int _dimension,
			const Norm_ptr_t &_norm,
			const constructDualityMapping_t &_constructDualityMapping);

	/** Getter for the this space.
	 *
	 * @return shared_ptr for this Space (correctly connected to the only
	 *         present instance via a weak_ptr), hence not as reference
	 */
	const NormedSpace_ptr_t getSpace() const
	{ return NormedSpace_ptr_t(Space); }

	/** Getter for the dual space.
	 *
	 * @return const reference to dual space
	 */
	const NormedSpace_ptr_t& getDualSpace() const
	{ return DualSpace; }

	/** Const getter for the norm in this space.
	 *
	 * @return norm object of this space
	 */
	const Norm_ptr_t& getNorm() const
	{ return norm; }

	/** Getter for the duality mapping
	 *
	 * @return const reference to duality mapping
	 */
	const Mapping_ptr_t& getDualityMapping() const;

	/** Const getter for the dimension of the representations in this space.
	 *
	 * @return dimension
	 */
	const unsigned int getDimension() const
	{ return dimension; }

	/** Creates a single element in this normed space.
	 *
	 * @param _dimension
	 * @return element instance
	 */
	SpaceElement_ptr_t createElement() const;

private:
	/** Weak reference to this space due to interface constraints.
	 * Otherwise this NormedSpace instance is never deleted as itself
	 * always holds a shared_ptr, keeping the reference count at 1.
	 */
	const boost::weak_ptr<NormedSpace> Space;

	//!> Reference to the dual of this space.
	const NormedSpace_ptr_t DualSpace;

	//!> Norm object
	const Norm_ptr_t norm;

	//!> stored function call for lazy creation of duality mapping
	const constructDualityMapping_t constructDualityMapping;

	//!> duality mapping from this into the dual space
	const Mapping_ptr_t dualitymapping;

	//!> dimension of the representations in this vector space
	const unsigned int dimension;

	//!> grant SpaceElement write access to opcounts
	friend class SpaceElement;

	//!> grant Norm write access to opcounts
	friend class Norm;

	//!> internal object for operator counting
	mutable VectorSpaceOperationCounts opcounts;

public:
	/** Getter for the operation counting structure.
	 *
	 * @return const ref to internal VectorSpaceOperationCounts
	 */
	const VectorSpaceOperationCounts& getOpCounts() const
	{ return opcounts; }
};


#endif /* NORMEDSPACE_HPP_ */
