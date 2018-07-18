/*
 * NormedDualSpace.hpp
 *
 *  Created on: Feb 5, 2015
 *      Author: heber
 */

#ifndef MINIMIZATIONS_SPACES_NORMEDDUALSPACE_HPP_
#define MINIMIZATIONS_SPACES_NORMEDDUALSPACE_HPP_

#include "BassoConfig.h"

#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/types.hpp"

/** This is a specialization of NormedSpace to break the shared_ptr cycle
 * between a space and its dual space. If both point at each other with
 * shared_ptr, then none will ever get deallocated. Hence, one of them
 * needs to point just via a weak_ptr. This is the DualSpace.
 */
class NormedDualSpace : public NormedSpace
{
private:
	//!> grant factory access to default cstor.
	friend class NormedSpaceFactory;

	/** Private default constructor to which only NormedSpaceFactory has access.
	 *
	 * @param _dimension dimension of this space
	 */
	NormedDualSpace(
			const unsigned int _dimension) :
		NormedSpace(_dimension)
	{}

	/** Setter of the dual space.
	 *
	 * This is for the factory only to allow it to connect the various
	 * spaces among one another.
	 *
	 * @param _dualspace
	 */
	void setDualSpace(const NormedSpace_ptr_t& _dualspace)
	{	const_cast<NormedSpace_weakptr_t&>(DualSpace) = _dualspace;	}

	/** Getter for the dual space.
	 *
	 * @return const reference to dual space
	 */
	const NormedSpace_ptr_t getDualSpace() const
	{ return NormedSpace_ptr_t(DualSpace); }

public:
	/** Constructor of NormedDualSpace.
	 *
	 * \warning setSpace() must be called to provide shared_ptr instace of
	 * this object. Otherwise getSpace() will return null and
	 * createElement() does not work, rendering this space unuseable.
	 *
	 * @param _dimension dimension of this space
	 * @param _norm norm object of the space
	 */
	NormedDualSpace(
			const unsigned int _dimension,
			const Norm_ptr_t &_norm) :
		NormedSpace(_dimension,_norm)
	{}

	virtual ~NormedDualSpace() {}

private:
	//!> here, we only have a weak_ptr to the dual space.
	const NormedSpace_weakptr_t DualSpace;
};



#endif /* MINIMIZATIONS_SPACES_NORMEDDUALSPACE_HPP_ */
