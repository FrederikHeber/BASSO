/*
 * Specifics/LInfinityDualityMapping.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef LINFINITYDUALITYMAPPING_HPP_
#define LINFINITYDUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"

class LInfinityDualityMappingUnitTest;

/** This class contains a duality mapping instance from a specific
 * l_infty space to its dual.
 *
 */
class LInfinityDualityMapping : public PowerTypeDualityMapping
{
	//!> grant unit test LInfinityDualityMappingUnitTest access to private parts.
	friend class LInfinityDualityMappingUnitTest;
public:
	/** Constructor for class LInfinityDualityMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _power power type of this duality mapping
	 */
	LInfinityDualityMapping(
			const NormedSpace_ptr_t &_NormedSpaceRef,
			const double _power) :
		PowerTypeDualityMapping(_NormedSpaceRef, _power),
		count(0),
		timing(boost::chrono::nanoseconds(0))
	{}

	/** Evaluates duality mapping at \a _x.
	 *
	 * \param _x point where to evaluate
	 */
	virtual const SpaceElement_ptr_t operator()(
			const SpaceElement_ptr_t &_x) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	const Mapping_ptr_t getAdjointMapping() const;

	/** Returns the number of times the operator() was called.
	 *
	 * @return number of calls
	 */
	const unsigned int getCount() const
	{ return count; }

	/** Returns the total runtime the program spent so far on
	 * its operator().
	 *
	 * @return runtime summed over all calls
	 */
	const boost::chrono::nanoseconds getTiming() const
	{ return timing; }

private:
	//!> number of times the operator was called
	mutable unsigned int count;

	//!> total runtime spent on executing this operator
	mutable boost::chrono::nanoseconds timing;
};

#endif /* LINFINITYDUALITYMAPPING_HPP_ */
