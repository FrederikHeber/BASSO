/*
 * Minimization_common.hpp
 *
 *  Created on: Nov 5, 2014
 *      Author: heber
 */

#ifndef MINIMIZATIONS_FUNCTIONS_MINIMIZATION_COMMON_HPP_
#define MINIMIZATIONS_FUNCTIONS_MINIMIZATION_COMMON_HPP_

#include "BassoConfig.h"

namespace Minimization
{

	//!> enumerates possible return codes of \sa CheckGradient().
	enum GradientStatus {
		gradient_success=0,	//!< gradient is (almost) zero
		gradient_continue=1,	//!< gradient is not yet zero
		error_noprogress=2, //!< gradient has not changed
		MAX_GradientStatus
	};

};



#endif /* MINIMIZATIONS_FUNCTIONS_MINIMIZATION_COMMON_HPP_ */
