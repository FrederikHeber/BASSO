/*
 * VectorProjection.hpp
 *
 *  Created on: Sep 4, 2014
 *      Author: heber
 */

#ifndef VECTORPROJECTION_HPP_
#define VECTORPROJECTION_HPP_

#include "BassoConfig.h"

#include <utility>

#include "Minimizations/types.hpp"

class PowerTypeDualityMapping;
class Norm;

/** This functor implements the projection of a vector onto another in
 * an (Lp-) Banach space.
 */
class VectorProjection
{
public:
	VectorProjection(
			const Norm &_lpnorm,
			const PowerTypeDualityMapping &_J_p,
			const double _p
			);
	~VectorProjection() {}

	/** Returns minimum, minimizer for the minimization functional depending
	 * on the Bregman distance to the scaled vector \a _projectonto.
	 *
	 * Basically, the minimizer is the length of \a _projectonto such that
	 * is the Bregman projection of \a _tobeprojected onto \a _projectonto.
	 *
	 * \warning Bregman distance is not symmetric. Hence, it matters which vector
	 * is projected onto which.
	 *
	 * @param _projectonto vector to project onto
	 * @param _tobeprojected vector to project
	 * @param _Tol tolerance threshold of minimization
	 * @return
	 */
	const std::pair<double, double> operator()(
			const SpaceElement_ptr_t &_projectonto,
			const SpaceElement_ptr_t &_tobeprojected,
			const double _Tol
			) const;

private:
	//!> lp Norm object
	const Norm &lpnorm;
	//!> DualityMapping object
	const PowerTypeDualityMapping &J_p;
	//!> value p of the Lp norm
	const double p;
};


#endif /* VECTORPROJECTION_HPP_ */
