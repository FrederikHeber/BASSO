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

SpaceElement::SpaceElement(const NormedSpace_weakptr_t _ref) :
	NormedSpaceRef(_ref),
	vector(Eigen::VectorXd::Zero(NormedSpace_ptr_t(_ref)->getDimension()))
{}

const bool SpaceElement::isInSpace(
		const NormedSpace_ptr_t &_NormedSpaceRef
		) const
{
	return _NormedSpaceRef == getSpace();
}

const bool SpaceElement::isZero(const double _threshold) const
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorComparison>(
					getSpace()->opcounts.instance));
	return vector.isZero(_threshold);
}

const bool SpaceElement::isApproxToConstant(
		const double _constant,
		const double _tolerance) const
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorComparison>(
					getSpace()->opcounts.instance));
	return vector.isApproxToConstant(_constant, _tolerance);
}

const bool SpaceElement::isApprox(
		const SpaceElement_ptr_t &_other,
		const double _tolerance) const
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorComparison>(
					getSpace()->opcounts.instance));
	return vector.isApprox(_other->vector, _tolerance);
}

const bool SpaceElement::isApprox(
		const SpaceElement &_other,
		const double _tolerance) const
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorComparison>(
					getSpace()->opcounts.instance));
	return vector.isApprox(_other.vector, _tolerance);
}

const bool SpaceElement::isNonnegative() const
{
	unsigned int rowMax;
	unsigned int colMax;
	const double value = vector.minCoeff(&rowMax, &colMax);
	return (value >= 0);
}

const SpaceElement_ptr_t SpaceElement::getSignVector() const
{
	SpaceElement_ptr_t signvector(getSpace()->createElement());
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorModification>(
					getSpace()->opcounts.instance));
	*signvector = Helpers::signum(vector);
	return signvector;
}

const SpaceElement_ptr_t SpaceElement::getAbsVector() const
{
	SpaceElement_ptr_t absvector(getSpace()->createElement());
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorModification>(
					getSpace()->opcounts.instance));
	*absvector = vector.array().abs();
	return absvector;
}

const SpaceElement_ptr_t SpaceElement::getCircShiftedVector(
		const int shift) const
{
	SpaceElement_ptr_t shiftedvector(getSpace()->createElement());
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorModification>(
					getSpace()->opcounts.instance));
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

const double& SpaceElement::operator[](const int i) const
{
	return vector[i];
}

const double SpaceElement::Norm() const
{
	return getSpace()->getNorm()->operator()(
			SpaceElement_ptr_t(SelfRef));
}

void SpaceElement::scaledAddition(
		const double _alpha,
		const SpaceElement_ptr_t &_element)
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorAddition>(
					getSpace()->opcounts.instance));
	assert( getSpace() == _element->getSpace() );
	assert( vector.innerSize() == _element->vector.innerSize() );
	vector += _alpha * _element->vector;
}

void SpaceElement::scaledAddition(
		const double _alpha,
		const SpaceElement &_element)
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorAddition>(
					getSpace()->opcounts.instance));
	assert( getSpace() == _element.getSpace() );
	assert( vector.innerSize() == _element.vector.innerSize() );
	vector += _alpha * _element.vector;
}

SpaceElement_ptr_t SpaceElement::operator*(const double _alpha) const
{
	SpaceElement_ptr_t newelement(getSpace()->createElement());
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::ScalarVectorMultiplication>(
					getSpace()->opcounts.instance));
	*newelement = vector;
	*newelement *= _alpha;
	return newelement;
}

const double SpaceElement::operator*(const SpaceElement_ptr_t &_element) const
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorMultiplication>(
					getSpace()->opcounts.instance));
	return vector.transpose() * _element->vector;
}

const double SpaceElement::operator*(const SpaceElement &_element) const
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorMultiplication>(
					getSpace()->opcounts.instance));
	return vector.transpose() * _element.vector;
}

SpaceElement_ptr_t SpaceElement::operator+(const SpaceElement_ptr_t &_element) const
{
	SpaceElement_ptr_t newelement(getSpace()->createElement());
	// exclude element creation time
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorAddition>(
					getSpace()->opcounts.instance));
	*newelement = vector;
	*newelement += _element;
	return newelement;
}

SpaceElement_ptr_t SpaceElement::operator-(const SpaceElement_ptr_t &_element) const
{
	SpaceElement_ptr_t newelement(getSpace()->createElement());
	// exclude element creation time
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorAddition>(
					getSpace()->opcounts.instance));
	*newelement = vector;
	*newelement -= _element;
	return newelement;
}

SpaceElement_ptr_t SpaceElement::operator*=(const double _alpha)
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::ScalarVectorMultiplication>(
					getSpace()->opcounts.instance));
	vector *= _alpha;
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator+=(const SpaceElement_ptr_t &_element)
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorAddition>(
					getSpace()->opcounts.instance));
	assert( getSpace() == _element->getSpace() );
	assert( vector.innerSize() == _element->vector.innerSize() );
	vector += _element->vector;
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator-=(const SpaceElement_ptr_t &_element)
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorAddition>(
					getSpace()->opcounts.instance));
	assert( getSpace() == _element->getSpace() );
	assert( vector.innerSize() == _element->vector.innerSize() );
	vector -= _element->vector;
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator=(const SpaceElement_ptr_t &_element)
{
	// check self-assignment
	if (this != _element.get()) {
		TIMEKEEPER(VectorSpaceOperations::getCountTiming<
				VectorSpaceOperations::VectorAssignment>(
						getSpace()->opcounts.instance));
		assert( getSpace() == _element->getSpace() );
		assert( vector.innerSize() == _element->vector.innerSize() );
		vector = _element->vector;
	}
	return SpaceElement_ptr_t(SelfRef);
}

SpaceElement_ptr_t SpaceElement::operator=(const Eigen::VectorXd &_vector)
{
	TIMEKEEPER(VectorSpaceOperations::getCountTiming<
			VectorSpaceOperations::VectorAssignment>(
					getSpace()->opcounts.instance));
	assert( vector.innerSize() == _vector.innerSize() );
	vector = _vector;
	return SpaceElement_ptr_t(SelfRef);
}
bool SpaceElement::operator==(const SpaceElement &_element)
{
	if (getSpace() != _element.getSpace())
		return false;
	if (vector != _element.vector)
		return false;
	return true;
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

