/*
 * BregmanDistance.hpp
 *
 *  Created on: Apr 25, 2014
 *      Author: heber
 */

#ifndef BREGMANDISTANCE_HPP_
#define BREGMANDISTANCE_HPP_

#include "BassoConfig.h"

#include "Minimizations/types.hpp"

class Mapping;
class Norm;

/** This implements a functor calculating the Bregman distance between
 * two points in a Lp space.
 *
 */
class BregmanDistance
{
public:
	/** Constructor for class BregmanDistance.
	 *
	 * @param _norm ref to norm
	 * @param _J_p ref to duality mapping
	 * @param _power power of duality mapping
	 */
	BregmanDistance(
			const Norm &_norm,
			const Mapping &_J_p,
			const double _power);

	/** Constructor for class BregmanDistance.
	 *
	 * @param _problem inverse problem with refs to norm and duality mapping
	 */
	BregmanDistance(
			const InverseProblem_ptr_t &_problem);

	~BregmanDistance() {}

	/** Calculate the Bregman distance between \a _x and \a _y.
	 *
	 * The Bregman distance reads as
	 * \f$ \Delta_p(x,y) = \frac 1 q ||x||^p + \frac 1 p ||y||^p - \langle J_{p}(x), y \rangle \f$
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
	 * @return Bregman distance between first and second argument
	 */
	double operator()(
			const SpaceElement_ptr_t &_x,
			const SpaceElement_ptr_t &_y
			) const;

	/** Calculate the Bregman distance between \a _x and \a _y with dual \a _xdual.
	 *
	 * The Bregman distance reads as
	 * \f$ \Delta_p(x,y) = \frac 1 q ||x||^p + \frac 1 p ||y||^p - \langle x^\ast, y \rangle \f$
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
	 * @param _xdual dual element to the first argument
	 * @return Bregman distance between first and second argument
	 */
	double operator()(
			const SpaceElement_ptr_t &_x,
			const SpaceElement_ptr_t &_y,
			const SpaceElement_ptr_t &_xdual
			) const;

private:
	//!> power type of the weight function of the duality mapping \a J_p
	const double power;
	//!> lp Norm object
	const Norm &norm;
	//!> LpDualityMapping object
	const Mapping &J_p;
};


#endif /* BREGMANDISTANCE_HPP_ */
