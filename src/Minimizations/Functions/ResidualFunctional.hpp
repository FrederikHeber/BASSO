/*
 * ResidualFunctional.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef RESIDUALFUNCTIONAL_HPP_
#define RESIDUALFUNCTIONAL_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/LpDualityMapping.hpp"
#include "Minimizations/Minimizers/LandweberMinimizer.hpp"

class ResidualFunctional
{
public:
	ResidualFunctional(
			const InverseProblem_ptr_t &_problem,
			const SpaceElement_ptr_t &_dualx,
			const SpaceElement_ptr_t &_u,
			const GeneralMinimizer &_mininizer
			) :
		problem(_problem),
		dualx(_dualx),
		u(_u),
		mininizer(_mininizer),
		residualvector(_problem->y->getSpace()->createElement())
	{}

	double operator()(double _arg) const;

private:
	const InverseProblem_ptr_t &problem;
	const SpaceElement_ptr_t &dualx;
	const SpaceElement_ptr_t &u;
	const GeneralMinimizer &mininizer;
	mutable SpaceElement_ptr_t residualvector;
};



#endif /* RESIDUALFUNCTIONAL_HPP_ */
