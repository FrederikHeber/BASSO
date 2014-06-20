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
			const LpNorm &_lpnorm,
			const DualityMapping &_J_p,
			const double _p) :
				p(_p),
				lpnorm(_lpnorm),
				J_p(_J_p)
	{}
	~BregmanDistance() {}

	/** Calculate the Bregman distance between \a _x and \a _y.
	 *
	 * The Bregman distance reads as
	 * \f$ \Delta_p(x,y) = \tfrac 1 q ||x||^p + \tfrac 1 p ||y||^p - \langle J_{p}(x), y \rangle \f$
	 *
	 * Note that the argument \a p (and \a q with it) is ambigious because
	 * the lp-norm also has an argument p. For clarification, the p and q
	 * in the above formulas refer to the power type \a p of the duality
	 * mapping \f$ L_p \f$. Hence, q is conjugate to \a _power. We refer
	 * to  \f$ L_p \f$'s always as \a _power to distinguish clearly from
	 * the norm.
	 *
	 * That's also why \a p is set in the cstor where we create the norm
	 * object while the power type may be given in operator().
	 *
	 * @param _x first argument
	 * @param _y second argument
	 * @param _power power of duality mapping
	 * @return Bregman distance between first and second argument
	 */
	double operator()(
			const Eigen::VectorXd &_x,
			const Eigen::VectorXd &_y,
			const double _power
			) const;

private:
	//!> value p of the Lp norm
	const double p;
	//!> lp Norm object
	const LpNorm &lpnorm;
	//!> DualityMapping object
	const DualityMapping &J_p;
};


#endif /* BREGMANDISTANCE_HPP_ */
