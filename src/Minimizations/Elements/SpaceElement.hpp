/*
 * SpaceElement.hpp
 *
 *  Created on: Oct 23, 2014
 *      Author: heber
 */

#ifndef SPACEELEMENT_HPP_
#define SPACEELEMENT_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

#include <boost/weak_ptr.hpp>
#include <utility>
#include <Eigen/Dense>
#include <iosfwd>
#include <limits>

class ElementCreator;
class NormedSpace;
class RepresentationAdvocate;

/** This class defines an element in a normed Space.
 *
 * Especially, various algebraic manipulations can be done with this
 * element or its underlying vector representation, respectively.
 *
 */
class SpaceElement
{
	//!> allow advocate access to private representation
	friend class RepresentationAdvocate;

	//!> allow NormedSpace access to private cstor
	friend class NormedSpace;

	//!> allow ElementCreator access to private representation
	friend class ElementCreator;

	//!> grant output operator access to private representation
	friend std::ostream & operator<<(
			std::ostream &ost,
			const SpaceElement &_element
			);

	//!> grant vector multiplication access to private representation
	friend const double operator*(
			const SpaceElement_ptr_t &_element,
			const SpaceElement_ptr_t &_otherelement);


	/** Private cstor to allow only a specific Space to create its elements.
	 *
	 * @param _ref reference to the NormedSpace this element belongs
	 */
	SpaceElement(const NormedSpace_weakptr_t _ref);

	/** Setter for the internal weak_ptr to \a this.
	 *
	 * \note This must be used by the NormedSpace::createElement() as
	 * otherwise certain operators will return shared_ptr's containing NULL.
	 *
	 * @param _selfref reference to this entity wrapped in weak_ptr.
	 */
	void setSelfRef(const SpaceElement_ptr_t &_selfref)
	{ const_cast<boost::weak_ptr<SpaceElement> &>(SelfRef) = _selfref; }

public:

	/** Checks whether this element is in the given NormedSpace.
	 *
	 * @param _NormedSpaceRef space to check
	 * @return true - element is in the given space, false - not
	 */
	const bool isInSpace(
			const NormedSpace_ptr_t &_NormedSpaceRef
			) const;

	/** Checks whether this contained representation is equal to zero.
	 *
	 * @param _threshold - tolerance for zero
	 * @return true - vector is zero, else - not
	 */
	const bool isZero(
			const double _threshold = std::numeric_limits<double>::epsilon()
			) const;

	/** Checks whether the representation is component-wise approximately
	 * equal the given \a _constant.
	 *
	 * @param _constant constant to check against
	 * @param _tolerance tolerance value
	 * @return true - every component is within [c-tol, c+tol], else - not
	 */
	const bool isApproxToConstant(
			const double _constant,
			const double _tolerance) const;

	/** Checks whether the representation is component-wise approximately
	 * equal the given \a _constant.
	 *
	 * \note While isApproxToConstant() checks on absolute precision,
	 * here we use relative precision comparing to the minimum of
	 * either component.
	 *
	 * @param _other other vector to compare to
	 * @param _tolerance tolerance value
	 * @return true - every component is within [c-tol, c+tol], else - not
	 */
	const bool isApprox(
			const SpaceElement_ptr_t &_other,
			const double _tolerance) const;

	/** Checks whether the representation is component-wise approximately
	 * equal the given \a _constant.
	 *
	 * \note While isApproxToConstant() checks on absolute precision,
	 * here we use relative precision comparing to the minimum of
	 * either component.
	 *
	 * @param _other other vector to compare to
	 * @param _tolerance tolerance value
	 * @return true - every component is within [c-tol, c+tol], else - not
	 */
	const bool isApprox(
			const SpaceElement &_other,
			const double _tolerance) const;

	/** Returns whether each component is non-negative.
	 *
	 * @return true - vector is non-negative, false - else
	 */
	const bool isNonnegative() const;

	/** Returns the vector with each component +1/0/-1.
	 *
	 * @return sign vector of this element
	 */
	const SpaceElement_ptr_t getSignVector() const;

	/** Returns the vector with each component as absolute value.
	 *
	 * @return abs vector of this element
	 */
	const SpaceElement_ptr_t getAbsVector() const;

	/** Circularly shifts the components of the representation by \a shift.
	 *
	 * @param shift number of entries to shift
	 * @return shifted copy of the vector
	 */
	const SpaceElement_ptr_t getCircShiftedVector(const int shift) const;

	/** Returns the index of the component with greatest absolute value.
	 *
	 * @return index of maximal component
	 */
	const std::pair<double, int> getMaxCoefficientAndIndex() const;

	/** Componentwise access to element's representation.
	 *
	 * @param i index of component
	 * @return \a i th component of element's representation.
	 */
	double& operator[](const int i);

	/** Componentwise access to element's representation.
	 *
	 * @param i index of component
	 * @return \a i th component of element's representation.
	 */
	const double& operator[](const int i) const;

	/** Calculates the norm of this element in the given Space.
	 *
	 * @return norm of the element
	 */
	const double Norm() const;

	/** Scaled addition.
	 *
	 * @param _alpha scalar factor
	 * @param _element element to scale and add
	 * @return instance plus \a _element that is scaled by \a _alpha
	 */
	void scaledAddition(const double _alpha, const SpaceElement_ptr_t &_element);

	/** Scaled addition.
	 *
	 * @param _alpha scalar factor
	 * @param _element element to scale and add
	 * @return instance plus \a _element that is scaled by \a _alpha
	 */
	void scaledAddition(const double _alpha, const SpaceElement &_element);

	/** Scalar multiplication.
	 *
	 * @param _alpha scalar factor
	 * @return new element scaled by \a_ alpha
	 */
	SpaceElement_ptr_t operator*(const double _alpha) const;

	/** Scalar multiplication.
	 *
	 * @param _alpha scalar factorprotected:
	 * @return reference to this
	 */
	SpaceElement_ptr_t operator*=(const double _alpha);

	/** Scalar product.
	 *
	 * @param _element other element
	 * @return scalar product between this and \a _element vector
	 */
	const double operator*(const SpaceElement_ptr_t &_element) const;

	/** Scalar product.
	 *
	 * @param _element other element
	 * @return scalar product between this and \a _element vector
	 */
	const double operator*(const SpaceElement &_element) const;

	/** Element sum.
	 *
	 * @param _element other element
	 * @return new element being the sum of this and \a other
	 */
	SpaceElement_ptr_t operator+(const SpaceElement_ptr_t &_element) const;

	/** Element sum.
	 *
	 * @param _element other element
	 * @return reference to this
	 */
	SpaceElement_ptr_t operator+=(const SpaceElement_ptr_t &_element);

	/** Element difference.
	 *
	 * @param _element other element
	 * @return new element being the difference of this and \a other
	 */
	SpaceElement_ptr_t operator-(const SpaceElement_ptr_t &_element) const;

	/** Element difference.
	 *
	 * @param _element other element
	 * @return reference to this
	 */
	SpaceElement_ptr_t operator-=(const SpaceElement_ptr_t &_element);

	/** Assignment operator
	 *
	 * @param _element other element
	 * @return reference to this
	 */
	SpaceElement_ptr_t operator=(const SpaceElement_ptr_t &_element);

	/** Equality operator
	 *
	 * @param _element other element
	 * @return true - both have equal spaces and components, false - else
	 */
	bool operator==(const SpaceElement &_element);

	/** Const getter to the space which this elements belongs to.
	 *
	 * @return const ref to space
	 */
	const NormedSpace_ptr_t getSpace() const
	{ return NormedSpace_ptr_t(NormedSpaceRef); }

	/** Sets the whole representation vector to zero.
	 *
	 */
	void setZero()
	{ vector.setZero(); }

private:
	/** Assignment operator for an Eigen::VectorXd.
	 *
	 * @param _vector new content of \a vector
	 * @return reference to this
	 */
	SpaceElement_ptr_t operator=(const Eigen::VectorXd &_vector);

	/** Const getter for the underlying vector representation.
	 *
	 * @return const ref to the vector representation
	 */
	const Eigen::VectorXd& getVectorRepresentation() const
	{ return vector; }

private:
	/** Similar as to internal ref for NormedSpace, the SpaceElements
	 * also holds a weak_ptr reference to itself in order to return
	 * it on specific operators
	 */
	const SpaceElement_weakptr_t SelfRef;

	//!> reference to space to which this element belongs
	const NormedSpace_weakptr_t NormedSpaceRef;

	//!> Vector object representation of the element
	Eigen::VectorXd vector;
};

/** Output iterator for this instance.
 *
 * @param ost stream to print to
 * @param _element element in space
 * @return ref to stream for concatenation
 */
std::ostream & operator<<(
		std::ostream &ost,
		const SpaceElement &_element
		);

/** Output iterator for this instance.
 *
 * @param ost stream to print to
 * @param _element element in space
 * @return ref to stream for concatenation
 */
std::ostream & operator<<(
		std::ostream &ost,
		const SpaceElement_ptr_t &_element
		);

/** Scalar multiplication.
 *
 * @param _alpha scalar factor
 * @param _element element in space
 * @return new element scaled by \a_ alpha
 */
inline SpaceElement_ptr_t operator*(
		const SpaceElement_ptr_t &_element,
		const double _alpha)
{
	return (*_element * _alpha);
}

/** Scalar multiplication.
 *
 * @param _element element in space
 * @param _alpha scalar factor
 * @return new element scaled by \a_ alpha
 */
inline SpaceElement_ptr_t operator*(
		const double _alpha,
		const SpaceElement_ptr_t &_element)
{
	return (*_element * _alpha);
}

/** Scalar product.
 *
 * @param _element element in space
 * @param _otherelement other element
 * @return scalar product between \a _element and \a _otherelement vector
 */
const double operator*(
		const SpaceElement_ptr_t &_element,
		const SpaceElement_ptr_t &_otherelement);

/** Element sum.
 *
 * @param _element first element
 * @param _otherelement other element
 * @return new element being the sum of this and \a other
 */
inline SpaceElement_ptr_t operator+(
		const SpaceElement_ptr_t &_element,
		const SpaceElement_ptr_t &_otherelement)
{
	return (*_element + _otherelement);
}

/** Element difference.
 *
 * @param _element first element
 * @param _otherelement other element
 * @return new element being the difference of this and \a other
 */
inline SpaceElement_ptr_t operator-(
		const SpaceElement_ptr_t &_element,
		const SpaceElement_ptr_t &_otherelement)
{
	return (*_element - _otherelement);
}

#endif /* SPACEELEMENT_HPP_ */
