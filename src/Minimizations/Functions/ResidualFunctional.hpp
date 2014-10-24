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

#include "Minimizations/LandweberMinimizer.hpp"
#include "Minimizations/DualityMappings/DualityMapping.hpp"

class ResidualFunctional
{
public:
	ResidualFunctional(
			const Eigen::VectorXd &_x,
			const Eigen::VectorXd &_dualx,
			const Eigen::VectorXd &_u,
			const Eigen::MatrixXd &_A,
			const Eigen::VectorXd &_y,
			const LandweberMinimizer &_landweber
			) :
				x(_x),
				dualx(_dualx),
				u(_u),
				A(_A),
				y(_y),
				landweber(_landweber)
	{}

	double operator()(double _arg) const
	{
		Eigen::VectorXd dual_solution = dualx;
		dual_solution -= _arg * u;
		Eigen::VectorXd x =
				landweber.J_q(dual_solution, landweber.DualPowerX);
		// calculate residual at candidate position (goes into hx[0])
		Eigen::VectorXd residual;
		return landweber.calculateResidual( x, A, y, residual );
	}

private:
	const Eigen::VectorXd &x;
	const Eigen::VectorXd &dualx;
	const Eigen::VectorXd &u;
	const Eigen::MatrixXd &A;
	const Eigen::VectorXd &y;
	const LandweberMinimizer &landweber;
};



#endif /* RESIDUALFUNCTIONAL_HPP_ */
