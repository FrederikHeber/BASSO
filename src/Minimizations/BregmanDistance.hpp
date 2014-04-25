/*
 * BregmanDistance.hpp
 *
 *  Created on: Apr 25, 2014
 *      Author: heber
 */

#ifndef BREGMANDISTANCE_HPP_
#define BREGMANDISTANCE_HPP_

#include "BassoConfig.h"

#include <Eigen/Dense>

#include "MinimizationExceptions.hpp"
#include "Minimizations/DualityMapping.hpp"
#include "Minimizations/LpNorm.hpp"

/** This implements a functor calculating the Bregman distance between
 * two points in a Lp space.
 *
 */
class BregmanDistance
{
public:
	BregmanDistance(
			const double _p) :
				p(_p),
				q(p/(p-1.)),
				lpnorm(p),
				J_p(p)
	{
		if (p <= 1.)
			throw MinimizationIllegalValue_exception()
				<< MinimizationIllegalValue_name("p");
	}
	~BregmanDistance() {}

	/** Calculate the Bregman distance between \a _x and \a _y.
	 *
	 * The Bregman distance reads as
	 * \f$ \Delta_n(x,y) = \tfrac 1 q ||x||^p + \tfrac 1 p ||y||^p - \langle J_{\mathrm{power}}(x), y \rangle \f$
	 *
	 * @param _x first argument
	 * @param _y second argument
	 * @return Bregman distance between first and second argument
	 */
	double operator()(
			const Eigen::VectorXd &_x,
			const Eigen::VectorXd &_y
			) const;

private:
	//!> value p of the Lp norm
	const double p;
	//!> value q that is conjugate to p
	const double q;
	//!> lp Norm object
	LpNorm lpnorm;
	//!> DualityMapping object
	DualityMapping J_p;
};


#endif /* BREGMANDISTANCE_HPP_ */
