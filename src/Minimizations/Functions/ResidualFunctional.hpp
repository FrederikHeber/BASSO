/*
 * ResidualFunctional.hpp
 *
 *  Created on: Oct 24, 2014
 *      Author: heber
 */

#ifndef RESIDUALFUNCTIONAL_HPP_
#define RESIDUALFUNCTIONAL_HPP_

#include "BassoConfig.h"

#include <boost/function.hpp>
#include <Eigen/Dense>

#include "Minimizations/InverseProblems/InverseProblem.hpp"
#include "Minimizations/Mappings/Specifics/LpDualityMapping.hpp"

class ResidualFunctional
{
public:
	//!> typedef for the residual calculation function
	typedef boost::function<double (
			const InverseProblem_ptr_t&,
			SpaceElement_ptr_t &)> calculateResidual_t;

	ResidualFunctional(
			const InverseProblem_ptr_t &_problem,
			const SpaceElement_ptr_t &_dualx,
			const SpaceElement_ptr_t &_u,
			const calculateResidual_t &_residualizer
			) :
		problem(_problem),
		dualx(_dualx),
		u(_u),
		residualizer(_residualizer),
		residualvector(_problem->y->getSpace()->createElement()),
		dual_solution(dualx->getSpace()->createElement())
	{}

	double operator()(double _arg) const;

private:
	const InverseProblem_ptr_t &problem;
	const SpaceElement_ptr_t &dualx;
	const SpaceElement_ptr_t &u;
	const calculateResidual_t &residualizer;
	mutable SpaceElement_ptr_t residualvector;
	mutable SpaceElement_ptr_t dual_solution;
};



#endif /* RESIDUALFUNCTIONAL_HPP_ */
