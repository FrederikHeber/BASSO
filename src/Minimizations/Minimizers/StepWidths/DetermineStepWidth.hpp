/*
 * DetermineStepWidth.hpp
 *
 *  Created on: Oct 30, 2014
 *      Author: heber
 */

#ifndef DETERMINESTEPWIDTH_HPP_
#define DETERMINESTEPWIDTH_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

/** Interface for step width calculation.
 *
 * Instantiate via DetermineStepWidthFactory.
 *
 */
struct DetermineStepWidth
{
	virtual const double operator()(
			const SpaceElement_ptr_t &_dualx,
			const SpaceElement_ptr_t &_u,
			const SpaceElement_ptr_t &_solution,
			const SpaceElement_ptr_t &_residual,
			const double _residuum,
			const double _TolX,
			const double _alpha
			) const = 0;
};


#endif /* DETERMINESTEPWIDTH_HPP_ */
