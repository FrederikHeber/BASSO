/*
 * RegularizedL1DualityMapping.hpp
 *
 *  Created on: May 01, 2017
 *      Author: heber
 */

#ifndef REGULARIZEDL1DUALITYMAPPING_HPP_
#define REGULARIZEDL1DUALITYMAPPING_HPP_

#include "BassoConfig.h"

#include <boost/chrono.hpp>

#include "Minimizations/Mappings/Specifics/L1DualityMapping.hpp"

class RegularizedL1DualityMappingUnitTest;

/** This class contains a duality mapping instance from a specific
 * l_1 space to its dual for the norm
 * \f$ \sqrt{ ||x||^2_1 + \lambda ||x||^2_2 } \f$.
 *
 */
class RegularizedL1DualityMapping : public L1DualityMapping
{
	//!> grant unit test RegularizedL1DualityMappingUnitTest access to private parts.
	friend class RegularizedL1DualityMappingUnitTest;
public:

	/** Constructor for class RegularizedL1DualityMapping.
	 *
	 * @param _NormedSpaceRef reference to space
	 * @param _power power type of this duality mapping
	 * @param _lambda regularization parameter
	 */
	RegularizedL1DualityMapping(
			const NormedSpace_weakptr_t &_NormedSpaceRef,
			const double _lambda = 0.1) :
		L1DualityMapping(_NormedSpaceRef, 2.),
		lambda(_lambda)
	{}

	//!> expose overloaded operator method from base class
	using Mapping::operator();

	/** Evaluates duality mapping at \a _x.
	 *
	 * \param _x point where to evaluate
	 * \param _Jx duality mapped \a _x
	 */
	virtual void operator()(
			const SpaceElement_ptr_t &_x,
			SpaceElement_ptr_t &_Jx) const;

	/** Gets the one dual element of possibly many that minimizes the
	 * scalar product to the given other element.
	 *
	 * \param _x point for which to get dual element(s)
	 * \param _y other element to minimize scalar product against
	 * \return _Jx duality mapped \a _x for which \f$ < Jx, y > \f$ is minimal
	 */
	virtual void getMinimumInfimum(
			const SpaceElement_ptr_t &_x,
			const SpaceElement_ptr_t &_y,
			SpaceElement_ptr_t &_Jx) const;

	/** Creates the adjoint mapping to this mapping.
	 *
	 * @return mapping instance with adjoint
	 */
	const Mapping_ptr_t getAdjointMapping() const;

	/** Setter for soft thresholding parameter \a lambda.
	 *
	 * @param _lambda new value for parameter
	 */
	void setLambda(const double _lambda) const
	{
		const_cast<double &>(lambda) = _lambda;
	}

	/** Getter for the regularization parameter.
	 *
	 * @return regularization parameter
	 */
	const double getLambda() const
	{ return lambda; }

protected:
	//!> soft thresholding parameter
	const double lambda;
};

#endif /* REGULARIZEDL1DUALITYMAPPING_HPP_ */
