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

#include "Math/Helpers.hpp"

#include "Minimizations/Norms/Norm.hpp"
#include "Minimizations/Spaces/NormedSpace.hpp"
#include "Minimizations/Spaces/VectorSpaceOperationCounts.hpp"

typedef VectorSpaceOperationCounts::TimeKeeper TimeKeeper;

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

const bool SpaceElement::isZero(const double _threshold) const
{
	TimeKeeper(NormedSpaceRef->opcounts.VectorComparison);
	return vector.isZero(_threshold);
}

const bool SpaceElement::isApproxToConstant(
		const double _constant,
		const double _tolerance) const
{
	TimeKeeper(NormedSpaceRef->opcounts.VectorComparison);
	return vector.isApproxToConstant(_constant, _tolerance);
}

const bool SpaceElement::isApprox(
		const SpaceElement_ptr_t &_other,
		const double _tolerance) const
{
	TimeKeeper(NormedSpaceRef->opcounts.VectorComparison);
	return vector.isApprox(_other->vector, _tolerance);
}

const bool SpaceElement::isApprox(
		const SpaceElement &_other,
		const double _tolerance) const
{
	TimeKeeper(NormedSpaceRef->opcounts.VectorComparison);
	return vector.isApprox(_other.vector, _tolerance);
}

const SpaceElement_ptr_t SpaceElement::getSignVector() const
{
	SpaceElement_ptr_t signvector(NormedSpaceRef->createElement());
	TimeKeeper(NormedSpaceRef->opcounts.VectorModification);
	*signvector = Helpers::signum(vector);
	return signvector;
}

const SpaceElement_ptr_t SpaceElement::getAbsVector() const
{
	SpaceElement_ptr_t absvector(NormedSpaceRef->createElement());
	TimeKeeper(NormedSpaceRef->opcounts.VectorModification);
	*absvector = vector.array().abs();
	return absvector;
}

const SpaceElement_ptr_t SpaceElement::getCircShiftedVector(
		const int shift) const
{
	SpaceElement_ptr_t shiftedvector(NormedSpaceRef->createElement());
	TimeKeeper(NormedSpaceRef->opcounts.VectorModification);
	*shiftedvector = Helpers::circshift(vector, shift);
	return shiftedvector;
}

const std::pair<double, int> SpaceElement::getMaxCoefficientAndIndex(
		) const
{
	unsigned int rowMax;
	unsigned int colMax;
	SpaceElement_ptr_t absvector = getAbsVector();
	const double value = absvector->vector.maxCoeff(&rowMax, &colMax);
	return std::make_pair(value, rowMax);
}

double& SpaceElement::operator[](const int i)
{
	return vector[i];
}

const double SpaceElement::operator[](const int i) const
{
	return vector[i];
}

const double SpaceElement::Norm() const
{
	return NormedSpaceRef->getNorm()->operator()(
			SpaceElement_ptr_t(SelfRef));
}

SpaceElement_ptr_t SpaceElement::operator*(const double _alpha) const
{
	SpaceElement_ptr_t newelement(NormedSpaceRef->createElement());
	TimeKeeper(NormedSpaceRef->opcounts.ScalarVectorMultiplication);
	*newelement = vector;
	*newelement *= _alpha;
	return newelement;
}

const double SpaceElement::operator*(const SpaceElement_ptr_t &_element) const
{
	TimeKeeper(NormedSpaceRef->opcounts.VectorMultiplication);
	return vector.transpose() * _element->vector;
}

const double SpaceElement::operator*(const SpaceElement &_element) const
{
	TimeKeeper(NormedSpaceRef->opcounts.VectorMultiplication);
	return vector.transpose() * _element.vector;
}

SpaceElement_ptr_t SpaceElement::operator+(const SpaceElement_ptr_t &_element) const
{
	SpaceElement_ptr_t newelement(NormedSpaceRef->createElement());
	// exclude element creation time
	TimeKeeper(NormedSpaceRef->opcounts.VectorAddition);
	*newelement = vector;
	*newelement += _element;
	return newelement;
}

SpaceElement_ptr_t SpaceElement::operator-(const SpaceElement_ptr_t &_element) const
{
	SpaceElement_ptr_t newelement(NormedSpaceRef->createElement());
	// exclude element creation time
	TimeKeeper(NormedSpaceRef->opcounts.VectorAddition);
	*newelement = vector;
	*newelement -= _element;
	return newelement;
}

SpaceElement_ptr_t SpaceElement::operator*=(const double _alpha)
{
	TimeKeeper(NormedSpaceRef->opcounts.ScalarVectorMultiplication);
	vector *= _alpha;
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator+=(const SpaceElement_ptr_t &_element)
{
	TimeKeeper(NormedSpaceRef->opcounts.VectorAddition);
	assert( NormedSpaceRef == _element->NormedSpaceRef );
	assert( vector.innerSize() == _element->vector.innerSize() );
	vector += _element->vector;
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator-=(const SpaceElement_ptr_t &_element)
{
	TimeKeeper(NormedSpaceRef->opcounts.VectorAddition);
	assert( NormedSpaceRef == _element->NormedSpaceRef );
	assert( vector.innerSize() == _element->vector.innerSize() );
	vector -= _element->vector;
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator=(const SpaceElement_ptr_t &_element)
{
	// check self-assignment
	if (this != _element.get()) {
		TimeKeeper(NormedSpaceRef->opcounts.VectorAssignment);
		assert( NormedSpaceRef == _element->NormedSpaceRef );
		assert( vector.innerSize() == _element->vector.innerSize() );
		vector = _element->vector;
	}
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator=(const Eigen::VectorXd &_vector)
{
	TimeKeeper(NormedSpaceRef->opcounts.VectorAssignment);
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
	return *_element * _otherelement;
}
