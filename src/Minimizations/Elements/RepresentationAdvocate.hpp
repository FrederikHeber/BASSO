/*
 * RepresentationAdvocate.hpp
 *
 *  Created on: Dec 11, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_ELEMENTS_REPRESENTATIONADVOCATE_HPP_
#define MINIMIZATIONS_ELEMENTS_REPRESENTATIONADVOCATE_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/Elements/SpaceElement.hpp"

class Backprojection;
class DiscretizedRadon;
class InverseProblemSolver;
class L1DualityMapping;
class LpDualityMapping;
struct LinearDependencyChecker;
class LinearMapping;
class LInfinityDualityMapping;
class RangeProjectionSolver;
class RelativeShrinkageMapping;
struct SingularValueDecomposition;
struct SpaceElementWriter;
struct SpaceElementReader;

/** The sole purpose of this class is to regulate the access to the
 * private representation of SpaceElement.
 *
 * I.e. herein we define all friends that have access to SpaceElement's
 * private representation.
 */
class RepresentationAdvocate
{
private:
	friend class Backprojection;
	friend class DiscretizedRadon;
	friend class InverseProblemSolver;
	friend class L1DualityMapping;
	friend class LpDualityMapping;
	friend struct LinearDependencyChecker;
	friend class LinearMapping;
	friend class LInfinityDualityMapping;
	friend class RangeProjectionSolver;
	friend class RelativeShrinkageMapping;
	friend struct SingularValueDecomposition;
	friend struct SpaceElementWriter;
	friend struct SpaceElementReader;
	template <class vector_type> friend struct VectorSetter;
	static const Eigen::VectorXd get(
			const SpaceElement_ptr_t & _element)
	{ return _element->getVectorRepresentation(); }

	static void set(
			SpaceElement_ptr_t & _element,
			const Eigen::VectorXd &_vector)
	{ 
		assert( _element->vector.innerSize() == _vector.innerSize() );
		_element->vector = _vector; 
	}
};


#endif /* MINIMIZATIONS_ELEMENTS_REPRESENTATIONADVOCATE_HPP_ */
