/*
 * HyperplaneProjection.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef HYPERPLANEPROJECTION_HPP_
#define HYPERPLANEPROJECTION_HPP_

#include "BassoConfig.h"

#include <vector>

#include "Minimizations/Functions/BregmanProjectionFunctional.hpp"
#include "Minimizations/Functions/Minimizers/MinimizationFunctional.hpp"
#include "Minimizations/types.hpp"

/** Structure containing all parameters to call BregmanProjectionFunctional functions.
 *
 * This is required to use function minimization that only allows
 * to pass a void* pointer to pass on information to the function to be
 * minimized.
 *
 * It calculates the projection onto the intersection of hyperplanes.
 *
 * \sa BregmanProjectionFunctional
 *
 */
struct HyperplaneProjection :
		public MinimizationFunctional< std::vector<double> >
{
	typedef typename MinimizationFunctional< std::vector<double> >::array_type array_type;

	BregmanProjectionFunctional &bregman;
	const SpaceElement_ptr_t &x;
	const std::vector<SpaceElement_ptr_t> &U;
	const std::vector<double> &alpha;

	/** Constructor to initialize refs.
	 *
	 */
	HyperplaneProjection(
		BregmanProjectionFunctional &_bregman,
		const SpaceElement_ptr_t &_x,
		const std::vector<SpaceElement_ptr_t> &_U,
		const std::vector<double> &_alpha
		);

	double function(const std::vector<double> &_value) const;

	const std::vector<double> gradient(const std::vector<double> &_value) const;

	void convertInternalTypeToArrayType(
			const std::vector<double> &_t,
			array_type & _x
			) const;

	void convertArrayTypeToInternalType(
			const array_type & _x,
			std::vector<double> &_t
			) const;
};



#endif /* HYPERPLANEPROJECTION_HPP_ */
