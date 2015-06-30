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
#include <limits>

/** Namespace Helpers contains helper functions for linear algebra routines
 *
 */
namespace Helpers {

	/** Created instance creates unique monotone increasing numbers.
	 *
	 */
	struct unique_number
	{
		unique_number() : number(0)
		{}

		unsigned int operator()()
		{ return number++; }

	private:
		unsigned int number;
	};


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

	/** This functor calculates the conjugate value q to a given p as
	 * \f$ q = p/(p-1.) \f$ such that \f$ 1/p + 1/q = 1 \f$.
	 *
	 * In case of infinity we always return 1 and the other way round.
	 *
	 * @param _p p value
	 * @return \f$ q = p/(p-1.) \f$ value
	 */
	inline double ConjugateValue(const double _p)
	{
		if (_p == 1.)
			return std::numeric_limits<double>::infinity();
		else if (_p == std::numeric_limits<double>::infinity())
			return  1.;
		else
			return _p/(_p-1.);
	}
};


#endif /* HELPERS_HPP_ */
