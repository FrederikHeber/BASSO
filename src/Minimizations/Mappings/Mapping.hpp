/*
 * Mapping.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef MAPPING_HPP_
#define MAPPING_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>

#include "Minimizations/types.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

/** Mapping represents a transformation of SpaceElements from one Space
 * into another.
 *
 * This basically defines the interface for all mappings.
 */
class Mapping
{
public:
	/** Constructor for a mapping.
	 *
	 * @param _SourceSpaceRef source space reference
	 * @param _TargetSpaceRef target space reference
	 */
	Mapping(
			const NormedSpace_weakptr_t &_SourceSpaceRef,
			const NormedSpace_weakptr_t &_TargetSpaceRef
			) :
		SourceSpaceRef(_SourceSpaceRef),
		TargetSpaceRef(_TargetSpaceRef)
	{}

	virtual ~Mapping() {}

	/** Getter for the source space of this mapping.
	 *
	 * @return ref to source space
	 */
	virtual const NormedSpace_ptr_t getSourceSpace() const
	{ return NormedSpace_ptr_t(SourceSpaceRef); }

	/** Getter for the target space of this mapping.
	 *
	 * @return ref to target space
	 */
	virtual const NormedSpace_ptr_t getTargetSpace() const
	{ return NormedSpace_ptr_t(TargetSpaceRef); }

	/** Dummy setter for tolerance.
	 *
	 * TODO: remove this as soon as PowerTypeDualityMapping is accessed directly again.
	 *
	 * \param _tolerance value to set to
	 */
	virtual void setTolerance(const double _tolerance) const
	{}

	/** Mapping function.
	 *
	 * @param _sourceelement element to map/transform
	 * @return new transformed/mapped element
	 */
	const SpaceElement_ptr_t operator()(
			const SpaceElement_ptr_t &_sourceelement
			) const
	{
		SpaceElement_ptr_t result = getTargetSpace()->createElement();
		operator()(_sourceelement, result);
		return result;
	}

	/** Mapping function.
	 *
	 * @param _sourceelement element to map/transform
	 * @param _destelement transformed/mapped element on return
	 */
	virtual void operator()(
			const SpaceElement_ptr_t &_sourceelement,
			SpaceElement_ptr_t &_destelement
			) const = 0;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	virtual const Mapping_ptr_t getAdjointMapping() const = 0;

	/** Returns the power type in case of a PowerTypeDualityMapping.
	 *
	 * @return 0 - not a PowerTypeDualityMapping, else - power type
	 */
	virtual const double getPower() const
	{ return 0.; }

	/** Returns the number of times the operator() was called.
	 *
	 * @return number of calls
	 */
	virtual const unsigned int getCount() const = 0;

	/** Returns the total runtime the program spent so far on
	 * its operator().
	 *
	 * @return runtime summed over all calls
	 */
	virtual const boost::chrono::nanoseconds getTiming() const = 0;

protected:
	//!> reference to space from which this mappings projects
	const NormedSpace_weakptr_t SourceSpaceRef;

	//!> reference to space into which this mappings projects
	const NormedSpace_weakptr_t TargetSpaceRef;
};



#endif /* MAPPING_HPP_ */
