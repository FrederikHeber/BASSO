/*
 * point_t.hpp
 *
 *  Created on: Jul 8, 2015
 *      Author: heber
 */

#ifndef COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_POINT_T_HPP_
#define COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_POINT_T_HPP_

#include "BassoConfig.h"

#include <iosfwd>

#include <Eigen/Dense>

/** Class for a point in the dim-dimensional space.
 *
 */
class point_t : public Eigen::Vector2d //boost::array<double, dim>
{
public:
	//!> constant for dimension of space
	enum { dim=2 };

	/** Default cstor for point_t.
	 *
	 */
	point_t();

	/** Conversion cstor from two doubles to point_t
	 *
	 * @param _x x component
	 * @param _y y component
	 */
	point_t(const double _x, const double _y);

	/** Conversion cstor from eigen vector to point_t
	 *
	 * @param _p eigen vector
	 */
	point_t(const Eigen::Vector2d &_p);

	/** Conversion cstor from eigen vector to point_t
	 *
	 * @param _p eigen vector
	 */
	point_t(const Eigen::Transpose<Eigen::Vector2d> &_p);

	/** Comparator for points in dim-dimensional space.
	 *
	 * @param _other other point
	 * @return true - point is lower than, to the left of \a _other
	 */
	bool operator<(const point_t &_other) const;
};

/** Output operator.
 *
 * @param _ost output stream to print to
 * @param _point point to print
 * @return ref to \a ost for concentation
 */
std::ostream& operator<<(
		std::ostream &_ost,
		const point_t &_point);

#endif /* COMPUTERTOMOGRAPHYBASE_DISCRETIZEDRADON_POINT_T_HPP_ */
