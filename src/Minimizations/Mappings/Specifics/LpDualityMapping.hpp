/*
 * Specifics/LpDualityMapping.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#ifndef LPDUALITYMAPPING_HPP_
#define LPDUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>

#include "Minimizations/Mappings/PowerTypeDualityMapping.hpp"
#include "Minimizations/Norms/LpNorm.hpp"

class LpDualityMappingUnitTest;

/** This class contains a duality mapping instance from a specific
 * l_p space with 1 < p < infty to its dual.
 *
 */
class LpDualityMapping : public PowerTypeDualityMapping
{
	//!> grant unit test LpDualityMappingUnitTest access to private parts.
	friend class LpDualityMappingUnitTest;
public:
	/** Constructor for class LpDualityMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _power power type of this duality mapping
	 */
	LpDualityMapping(
			const NormedSpace_weakptr_t &_NormedSpaceRef,
			const double _power);

	virtual ~LpDualityMapping() {}

	//!> expose overloaded operator method from base class
	using Mapping::operator();

	/** Mapping function.
	 *
	 * @param _sourceelement element to map/transform
	 * @param _destelement new transformed/mapped element
	 */
	void operator()(
			const SpaceElement_ptr_t &_sourceelement,
			SpaceElement_ptr_t &_destelement
			) const;

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

protected:
	//!> number of times the operator was called
	mutable unsigned int count;

	//!> total runtime spent on executing this operator
	mutable boost::chrono::nanoseconds timing;
};

#endif /* LPDUALITYMAPPING_HPP_ */
