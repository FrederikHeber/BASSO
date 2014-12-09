/*
 * SpaceElement.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: heber
 */


#include "BassoConfig.h"

#include "SpaceElement.hpp"

#include <cassert>
#include <iostream>

#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"

SpaceElement::SpaceElement(const NormedSpace_ptr_t &_ref) :
	NormedSpaceRef(_ref),
	vector(Eigen::VectorXd::Zero(_ref->getDimension()))
{}

const bool SpaceElement::isInSpace(
		const NormedSpace_ptr_t &_NormedSpaceRef
		) const
{
	return _NormedSpaceRef == getSpace();
}

const bool SpaceElement::isZero() const
{
	return vector.isZero();
}

const bool SpaceElement::isApproxToConstant(
		const double _constant,
		const double _tolerance) const
{
	return vector.isApproxToConstant(_constant, _tolerance);
}

const double SpaceElement::Norm() const
{
	return NormedSpaceRef->getNorm()->operator()(
			SpaceElement_ptr_t(SelfRef));
}

SpaceElement_ptr_t SpaceElement::operator*(const double _alpha) const
{
	SpaceElement_ptr_t newelement(NormedSpaceRef->createElement());
	*newelement = vector;
	*newelement *= _alpha;
	return newelement;
}

const double SpaceElement::operator*(const SpaceElement_ptr_t &_element) const
{
	return (vector.transpose() * _element->vector);
}

const double SpaceElement::operator*(const SpaceElement &_element) const
{
	return vector.transpose() * _element.vector;
}

SpaceElement_ptr_t SpaceElement::operator+(const SpaceElement_ptr_t &_element) const
{
	SpaceElement_ptr_t newelement(NormedSpaceRef->createElement());
	*newelement = vector;
	*newelement += _element;
	return newelement;
}

SpaceElement_ptr_t SpaceElement::operator-(const SpaceElement_ptr_t &_element) const
{
	SpaceElement_ptr_t newelement(NormedSpaceRef->createElement());
	*newelement = vector;
	*newelement -= _element;
	return newelement;
}

SpaceElement_ptr_t SpaceElement::operator*=(const double _alpha)
{
	vector *= _alpha;
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator+=(const SpaceElement_ptr_t &_element)
{
	assert( NormedSpaceRef == _element->NormedSpaceRef );
	assert( vector.innerSize() == _element->vector.innerSize() );
	vector += _element->vector;
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator-=(const SpaceElement_ptr_t &_element)
{
	assert( NormedSpaceRef == _element->NormedSpaceRef );
	assert( vector.innerSize() == _element->vector.innerSize() );
	vector -= _element->vector;
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator=(const SpaceElement_ptr_t &_element)
{
	// check self-assignment
	if (this != _element.get()) {
		assert( NormedSpaceRef == _element->NormedSpaceRef );
		assert( vector.innerSize() == _element->vector.innerSize() );
		vector = _element->vector;
	}
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator=(const Eigen::VectorXd &_vector)
{
	assert( vector.innerSize() == _vector.innerSize() );
	vector = _vector;
	return SpaceElement_ptr_t(SelfRef);
}

std::ostream & operator<<(std::ostream &ost, const SpaceElement &_element)
{
	ost << _element.getVectorRepresentation().transpose();
	return ost;
}

std::ostream & operator<<(std::ostream &ost, const SpaceElement_ptr_t &_element)
{
	ost << *_element;
	return ost;
}

/** Scalar product.
 *
 * @param _element element in space
 * @param _otherelement other element
 * @return scalar product between \a _element and \a _otherelement vector
 */
const double operator*(
		const SpaceElement_ptr_t &_element,
		const SpaceElement_ptr_t &_otherelement)
{
	return _element->getVectorRepresentation().transpose() *
			_otherelement->getVectorRepresentation();
}
