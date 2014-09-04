/*
 * Helpers.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#ifndef HELPERS_HPP_
#define HELPERS_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

/** Namespace Helpers contains helper functions for linear algebra routines
 *
 */
namespace Helpers {
	/** Implements a circular shift of the vector's components by \a shift.
	 *
	 * Note that the sign of \a shift decides on the direction.
	 *
	 * \param shift step-size (and direction) of the shift
	 * \return VectorXd with shifted components
	 */
	Eigen::VectorXd circshift(const Eigen::VectorXd &_x, const int _shift);

	/** Gives the sign (-1/1) of the value \a x.
	 *
	 * \param x value
	 * \return -1, 0, or 1 if \a x is negative, zero or positive
	 */
	double sign(const double _x);

	/** Implements signum function, i.e. componentwise sign.
	 *
	 * \param x vector
	 * \return x with -1, 0 or 1 per component
	 */
	Eigen::VectorXd signum(const Eigen::VectorXd &_x);
};


#endif /* HELPERS_HPP_ */
