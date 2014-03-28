/*
 * BregmanFunctional.hpp
 *
 *  Created on: Mar 27, 2014
 *      Author: heber
 */

#ifndef BREGMANFUNCTIONAL_HPP_
#define BREGMANFUNCTIONAL_HPP_

#include <Eigen/Dense>
#include <utility>

#include "Minimizations/DualityMapping.hpp"

/** Functor to calculate BregmanFunctional functional/distance.
 *
 */
template <unsigned int p>
class BregmanFunctional
{
public:
	BregmanFunctional(const double _tolerance) :
		tolerance(_tolerance)
	{}
	~BregmanFunctional() {}

	/** Implements BregmanFunctional functional.
	 *
	 * \param _t
	 * \param _x vector
	 * \param _U
	 * \param _alpha
	 * \param _q power of the weight of the duality mapping
	 */
	double operator()(
			const Eigen::VectorXd &_t,
			const Eigen::VectorXd &_x,
			const Eigen::MatrixXd &_U,
			const Eigen::VectorXd &_alpha,
			const unsigned int _q
			)
	{
		// x=x-U*t;
		const Eigen::VectorXd resx = _x - _U * _t;
		// fval=1/q*norm(x,p)^q+alpha'*t;
		const double fval =
				1./(double)_q * ::pow(resx.lpNorm<p>(), _q)
				+ _alpha.transpose() * _t;
		return fval;
	}

	/** Implements BregmanFunctional functional.
	 *
	 * \param _t
	 * \param _x vector
	 * \param _U
	 * \param _alpha
	 * \param _q power of the weight of the duality mapping
	 */
	Eigen::VectorXd gradient(
			const Eigen::VectorXd &_t,
			const Eigen::VectorXd &_x,
			const Eigen::MatrixXd &_U,
			const Eigen::VectorXd &_alpha,
			const unsigned int _q
			)
	{
		const double fval = operator()(_t,_x,_U_alpha,_q);
		const DualityMapping<p> J_p(_q);
		J_p.setTolerance(tolerance);
		const Eigen::VectorXd gval =
				_alpha -
				_U.transpose() * J_p(resx);

		return gval;
	}

private:
	//!> tolerance to pass on to used duality mapping instance
	const double tolerance;
};


#endif /* BREGMANFUNCTIONAL_HPP_ */
